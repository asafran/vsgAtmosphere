#version 450
#extension GL_ARB_separate_shader_objects : enable

#include "rendering_constants.glsl"
#include "constants.glsl"
#include "functions.glsl"

// ------------------------------------------------------------------
// CONSTANTS --------------------------------------------------------
// ------------------------------------------------------------------
layout (local_size_x_id = 1, local_size_y_id = 1, local_size_z = 1) in;
layout (constant_id = 2) const int CUBE_SIZE = 128;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 0) uniform sampler2D transmittance_texture;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 2) uniform sampler3D scattering_texture;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 3) uniform sampler3D single_mie_scattering_texture;

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 4, std140) uniform Settings
{
	vec4 whitePoint;
	vec2 sunSize;
} settings;

layout(set = POSITIONAL_DESCRIPTOR_SET, binding = 0, std140) uniform Positional
{
	vec3 sunDirection;
    float exposure;
    vec3 cameraPos;
    float radius;
    float mu_s;
} positional;

layout(set = PBR_DESCRIPTOR_SET, binding = 0, rgba16f) uniform imageCube prefilteredMap;


// ------------------------------------------------------------------
// FUNCTIONS --------------------------------------------------------
// ------------------------------------------------------------------

Luminance3 GetSolarRadiance()
{
  return solar_irradiance / (PI * sun_angular_radius * sun_angular_radius);
}

RadianceSpectrum GetSkyRadiance(IN(Direction) view_ray, IN(Number) nu, OUT(DimensionlessSpectrum) transmittance)
{
  // Compute the distance to the top atmosphere boundary along the view ray,
  // assuming the viewer is in space (or NaN if the view ray does not intersect
  // the atmosphere).
  Length r = positional.radius;
  Number rmu = dot(positional.cameraPos, view_ray);

  // Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
  Number mu = rmu / r;
  Number mu_s = positional.mu_s;
  bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

  transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum(0.0) :
      GetTransmittanceToTopAtmosphereBoundary(
           transmittance_texture, r, mu);
  IrradianceSpectrum single_mie_scattering;
  IrradianceSpectrum scattering = GetCombinedScattering(
         scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
  return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
      MiePhaseFunction(mie_phase_function_g, nu);
}
/*
// ------------------------------------------------------------------
// Based omn http://byteblacksmith.com/improvements-to-the-canonical-one-liner-glsl-rand-for-opengl-es-2-0/
float random(vec2 co)
{
	float a = 12.9898;
	float b = 78.233;
	float c = 43758.5453;
	float dt= dot(co.xy ,vec2(a,b));
	float sn= mod(dt,3.14);
	return fract(sin(sn) * c);
}

vec2 hammersley2d(uint i, uint N)
{
	// Radical inverse based on http://holger.dammertz.org/stuff/notes_HammersleyOnHemisphere.html
	uint bits = (i << 16u) | (i >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	float rdi = float(bits) * 2.3283064365386963e-10;
	return vec2(float(i) /float(N), rdi);
}

// Based on http://blog.selfshadow.com/publications/s2013-shading-course/karis/s2013_pbs_epic_slides.pdf
vec3 importanceSample_GGX(vec2 Xi, float roughness, vec3 normal)
{
	// Maps a 2D point to a hemisphere with spread based on roughness
	float alpha = roughness * roughness;
	float phi = 2.0 * PI * Xi.x + random(normal.xz) * 0.1;
	float cosTheta = sqrt((1.0 - Xi.y) / (1.0 + (alpha*alpha - 1.0) * Xi.y));
	float sinTheta = sqrt(1.0 - cosTheta * cosTheta);
	vec3 H = vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);

	// Tangent space
	vec3 up = abs(normal.z) < 0.999 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
	vec3 tangentX = normalize(cross(up, normal));
	vec3 tangentY = normalize(cross(normal, tangentX));

	// Convert to world Space
	return normalize(tangentX * H.x + tangentY * H.y + normal * H.z);
}

// Normal Distribution function
float D_GGX(float dotNH, float roughness)
{
	float alpha = roughness * roughness;
	float alpha2 = alpha * alpha;
	float denom = dotNH * dotNH * (alpha2 - 1.0) + 1.0;
	return (alpha2)/(PI * denom*denom);
}

vec3 prefilterEnvMap(vec3 R, float roughness)
{
	vec3 N = R;
	vec3 V = R;
	vec3 color = vec3(0.0);
	float totalWeight = 0.0;
	float envMapDim = float(textureSize(samplerEnv, 0).s);
	for(uint i = 0u; i < consts.numSamples; i++) {
		vec2 Xi = hammersley2d(i, consts.numSamples);
		vec3 H = importanceSample_GGX(Xi, roughness, N);
		vec3 L = 2.0 * dot(V, H) * H - V;
		float dotNL = clamp(dot(N, L), 0.0, 1.0);
		if(dotNL > 0.0) {
			// Filtering based on https://placeholderart.wordpress.com/2015/07/28/implementation-notes-runtime-environment-map-filtering-for-image-based-lighting/

			float dotNH = clamp(dot(N, H), 0.0, 1.0);
			float dotVH = clamp(dot(V, H), 0.0, 1.0);

			// Probability Distribution Function
			float pdf = D_GGX(dotNH, roughness) * dotNH / (4.0 * dotVH) + 0.0001;
			// Slid angle of current smple
			float omegaS = 1.0 / (float(consts.numSamples) * pdf);
			// Solid angle of 1 pixel across all cube faces
			float omegaP = 4.0 * PI / (6.0 * envMapDim * envMapDim);
			// Biased (+1.0) mip level for better result
			float mipLevel = roughness == 0.0 ? 0.0 : max(0.5 * log2(omegaS / omegaP) + 1.0, 0.0f);
			color += textureLod(samplerEnv, L, mipLevel).rgb * dotNL;
			totalWeight += dotNL;

		}
	}
	return (color / totalWeight);
}
*/
// ------------------------------------------------------------------

void main()
{
	vec2 coord = vec2(gl_GlobalInvocationID.xy) / float(CUBE_SIZE);
	uint layer = gl_GlobalInvocationID.z;

	vec3 view_direction = vec3(0.0, 0.0, 0.0);

	switch (layer) {
	case 0:
		view_direction = normalize(vec3(0.5, 0.5 - coord.y, coord.x - 0.5));
		break;
	case 1:
		view_direction = normalize(vec3(-0.5, 0.5 - coord.y, 0.5 - coord.x));
		break;
	case 2:
		view_direction = normalize(vec3(coord.x - 0.5, 0.5, 0.5 - coord.y));
		break;
	case 3:
		view_direction = normalize(vec3(coord.x - 0.5, -0.5, coord.y - 0.5));
		break;
	case 4:
		view_direction = normalize(vec3(coord.x - 0.5, 0.5 - coord.y, -0.5));
		break;
	case 5:
		view_direction = normalize(vec3(0.5 - coord.x, 0.5 - coord.y, 0.5));
		break;
	}

    Number nu = dot(view_direction, positional.sunDirection);

    // Compute the radiance of the sky.
    vec3 solarRadiance = GetSolarRadiance();
	vec3 transmittance;
	vec3 radiance = GetSkyRadiance(view_direction, nu, transmittance);
#ifdef ATMOSHPERE_RADIANCE
    radiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    solarRadiance *= SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
#endif

	// If the view ray intersects the Sun, add the Sun radiance.
	if (nu > settings.sunSize.y)
		radiance = radiance + transmittance * solarRadiance;

	ivec3 texel = ivec3(gl_GlobalInvocationID);
	imageStore(prefilteredMap, texel, vec4(pow(vec3(1.0) - exp(-radiance / settings.whitePoint.rbg * positional.exposure), vec3(1.0 / 2.2)), 1.0));
}


/*
void main()
{		
	vec3 N = normalize(inPos);
	outColor = vec4(prefilterEnvMap(N, consts.roughness), 1.0);
}
*/
