#version 450
#extension GL_ARB_separate_shader_objects : enable
#pragma import_defines (ATMOSHPERE_CLOUDS, ATMOSHPERE_VIEWER_IN_SPACE)

#include "rendering_constants.glsl"
#include "constants.glsl"
#include "functions.glsl"


layout(location = 0) in vec3 inRay;
layout(location = 1) in vec2 inUV;
layout(location = 0) out vec4 outColor;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 0) uniform sampler2D s_Transmittance;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 1) uniform sampler2D s_Irradiance;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 2) uniform sampler3D s_Scattering;
layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 3) uniform sampler3D s_SingleMieScattering;

layout(set = ATMOSHPERE_DESCRIPTOR_SET, binding = 4, std140) uniform Settings
{
	vec4 whitePoint;
	vec2 sunSize;
} settings;

layout(set = POSITIONAL_DESCRIPTOR_SET, binding = 0, std140) uniform Positional
{
	vec4 sunDirectionExp;
    vec4 cameraPosR;
    float mu_s;
} positional;

#ifdef ATMOSHPERE_CLOUDS
layout(set = CLOUDS_DESCRIPTOR_SET, binding = 0) uniform sampler3D s_ShapeNoise;
layout(set = CLOUDS_DESCRIPTOR_SET, binding = 1) uniform sampler3D s_DetailNoise;
layout(set = CLOUDS_DESCRIPTOR_SET, binding = 2) uniform sampler2D s_BlueNoise;
layout(set = CLOUDS_DESCRIPTOR_SET, binding = 3) uniform sampler2D s_CurlNoise;

const int NUM_CONE_SAMPLES = 6;

layout(set = CLOUDS_DESCRIPTOR_SET, binding = 4, std140) uniform Clouds
{
    float 	  cloudMinHeight;
    float 	  cloudMaxHeight;
    float 	  shapeNoiseScale;
    float 	  detailNoiseScale;

    float 	  detailNoiseModifier;
    float 	  turbulenceNoiseScale;
    float 	  turbulenceAmount;
    float 	  cloudCoverage;

    vec3 	  windDirection;
    float	  windSpeed;

    float	  windShearOffset;
    int 	  maxNumSteps;
    float 	  lightStepLength;
    float 	  lightConeRadius;

    vec3      cloudBaseColor;
    float 	  precipitation;

    float 	  ambientLightFactor;
    float 	  sunLightFactor;
    float 	  henyeyGreensteinGForward;
    float 	  henyeyGreensteinGBackward;
} clouds;
#endif

layout(push_constant) uniform PushConstants {
    mat4 projection;
    mat4 modelView;
} pc; //inverse

// ------------------------------------------------------------------
// FUNCTIONS --------------------------------------------------------
// ------------------------------------------------------------------


Luminance3 GetSolarRadiance() {
      return solar_irradiance /
          (PI * sun_angular_radius * sun_angular_radius) *
          SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
Luminance3 GetSkyRadiance(
        Position camera, Direction view_ray, Length shadow_length,
        Direction sun_direction, out DimensionlessSpectrum transmittance) {
      return GetSkyRadiance(s_Transmittance,
          s_Scattering, s_SingleMieScattering,
          camera, view_ray, shadow_length, sun_direction, transmittance) *
          SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
    Luminance3 GetSkyRadianceToPoint(
        Position camera, Position point, Length shadow_length,
        Direction sun_direction, out DimensionlessSpectrum transmittance) {
      return GetSkyRadianceToPoint(s_Transmittance,
          s_Scattering, s_SingleMieScattering,
          camera, point, shadow_length, sun_direction, transmittance) *
          SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
    Illuminance3 GetSunAndSkyIrradiance(
       Position p, Direction normal, Direction sun_direction,
       out IrradianceSpectrum sky_irradiance) {
      IrradianceSpectrum sun_irradiance = GetSunAndSkyIrradiance(
          s_Transmittance, s_Irradiance, p, normal,
          sun_direction, sky_irradiance);
      sky_irradiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
      return sun_irradiance * SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
#ifdef ATMOSHPERE_CLOUDS
// ------------------------------------------------------------------

struct Ray
{
    vec3 origin;
    vec3 direction;
};

// ------------------------------------------------------------------
// CLOUDS FUNCTIONS --------------------------------------------------------
// ------------------------------------------------------------------

Ray generate_ray()
{
    vec2 tex_coord_neg_to_pos = inUV * 2.0f - 1.0f;
    vec4 target = pc.projection * vec4(tex_coord_neg_to_pos, 0.0f, 1.0f);
    target /= target.w;

    Ray ray;

    ray.origin 	  = positional.cameraPosR.xyz;
    ray.direction = normalize(target.xyz - ray.origin);

    return ray;
}

// ------------------------------------------------------------------

vec3 ray_sphere_intersection(in Ray ray, in vec3 sphere_center, in float sphere_radius)
{
    vec3 l = ray.origin - sphere_center;
	float a = 1.0;
	float b = 2.0 * dot(ray.direction, l);
	float c = dot(l, l) - pow(sphere_radius, 2);
	float D = pow(b, 2) - 4.0 * a * c;

    if (D < 0.0)
		return ray.origin;
	else if (abs(D) - 0.00005 <= 0.0)
		return ray.origin + ray.direction * (-0.5 * b / a);
	else
	{
		float q = 0.0;
		if (b > 0.0)
            q = -0.5 * (b + sqrt(D));
		else
            q = -0.5 * (b - sqrt(D));

		float h1 = q / a;
		float h2 = c / q;
		vec2 t = vec2(min(h1, h2), max(h1, h2));

		if (t.x < 0.0)
        {
			t.x = t.y;

            if (t.x < 0.0)
				return ray.origin;
		}

		return ray.origin + t.x * ray.direction;
	}
}

// ------------------------------------------------------------------

float blue_noise()
{
	ivec2 size = textureSize(s_BlueNoise, 0);

	vec2 interleaved_pos = (mod(floor(gl_FragCoord.xy), float(size.x)));
	vec2 tex_coord 	     = interleaved_pos / float(size.x) + vec2(0.5f / float(size.x), 0.5f / float(size.x));

	return texture(s_BlueNoise, tex_coord).r * 2.0f - 1.0f;
}

// ------------------------------------------------------------------

float remap(float original_value, float original_min, float original_max, float new_min, float new_max)
{
	return new_min + (((original_value - original_min) / (original_max - original_min)) * (new_max - new_min));
}

// ------------------------------------------------------------------

// returns height fraction [0, 1] for point in cloud
float height_fraction_for_point(vec3 _position)
{
    float radius = length(_position);
    float height = (radius - bottom_radius ) / (top_radius - bottom_radius);
	return clamp(height, 0.0f, 1.0f);
}

// ------------------------------------------------------------------

float density_height_gradient_for_point(vec3 _position, float _height_fraction)
{
	return 1.0f;
}

// ------------------------------------------------------------------

float sample_cloud_density(vec3 _position, float _height_fraction, float _lod, bool _use_detail)
{
    // Shear cloud top along wind direction.
    vec3 position = _position + clouds.windDirection * clouds.windShearOffset * _height_fraction;

    // Animate clouds in wind direction and add a small upward bias to the wind direction.
    position += (clouds.windDirection + vec3(0.0f, 0.1f, 0.0f)) * clouds.windSpeed * 0.0; ///time

    // Read the low-frequency Perlin-Worley and Worley noises.
    vec4 low_frequency_noises = textureLod(s_ShapeNoise, position *  clouds.shapeNoiseScale, _lod);

    // Build an FBM out of the low-frequency Worley noises to add detail to the low-frequeny Perlin-Worley noise.
    float low_freq_fbm = (low_frequency_noises.g * 0.625f) + (low_frequency_noises.b * 0.25f) + (low_frequency_noises.a * 0.125f);

    // Define the base cloud shape by dilating it with the low-frequency FBM made of Worley noise.
    float base_cloud = remap(low_frequency_noises.r, (1.0f - low_freq_fbm), 1.0f, 0.0f, 1.0f);

    // Get the density-height gradient using the density height function.
    float density_height_gradient = density_height_gradient_for_point(position, _height_fraction);

    // Apply the height function to the base cloud shape.
    base_cloud *= density_height_gradient;

    // Remap to apply the cloud coverage attribute.
    float base_cloud_with_coverage = remap(base_cloud, clouds.cloudCoverage, 1.0f, 0.0f, 1.0f);

    // Multiply result by the cloud coverage attribute so that smaller clouds are lighter and more aesthetically pleasing.
    base_cloud_with_coverage *= clouds.cloudCoverage;

    // Exit out if base cloud density is zero.
    if (base_cloud_with_coverage <= 0.0f)
        return 0.0f;

    float final_cloud = base_cloud_with_coverage;

    if (_use_detail)
    {
        // Sample curl noise texture.
        vec2 curl_noise = textureLod(s_CurlNoise, position.xz * clouds.turbulenceNoiseScale, 0.0f).rg;

        // Add some turbulence to bottom of clouds.
        position.xy += curl_noise * (1.0f - _height_fraction) * clouds.turbulenceAmount;

        // Sample high-frequency noises.
        vec3 high_frequency_noises = textureLod(s_DetailNoise, position * clouds.detailNoiseScale, _lod).rgb;

        // Build high-frequency Worley noise FBM.
        float high_freq_fbm = (high_frequency_noises.r * 0.625f) + (high_frequency_noises.g * 0.25f) + (high_frequency_noises.b * 0.125f);

        // Transition from wispy shapes to billowy shapes over height.
        float high_freq_noise_modifier = mix(1.0f - high_freq_fbm, high_freq_fbm, clamp(_height_fraction * 10.0f, 0.0f, 1.0f));

        // Erode the base cloud shape with the distorted high-frequency Worley noise.
        final_cloud = remap(base_cloud_with_coverage, high_freq_noise_modifier * clouds.detailNoiseModifier, 1.0f, 0.0f, 1.0f);
    }

    return clamp(final_cloud, 0.0f, 1.0f);
}

// ------------------------------------------------------------------

float sample_cloud_density_along_cone(vec3 _position, vec3 _light_dir)
{
	const vec3 noise_kernel[6] =
	{
		{ -0.6, -0.8, -0.2 },
		{ 1.0, -0.3, 0.0 },
		{ -0.7, 0.0, 0.7 },
		{ -0.2, 0.6, -0.8 },
		{ 0.4, 0.3, 0.9 },
		{ -0.2, 0.6, -0.8 }
	};

	float density_along_cone = 0.0f;

	for (int i = 0; i < NUM_CONE_SAMPLES; i++)
	{
        // March ray forward along light direction.
		_position += _light_dir * clouds.lightStepLength;

        // Compute offset within the cone.
		vec3 random_offset = noise_kernel[i] * clouds.lightStepLength * clouds.lightConeRadius * (float(i + 1));

        // Add offset to position.
		vec3 p = _position + random_offset;

        // Compute height fraction for the current position.
		float height_fraction = height_fraction_for_point(p);

		// Skipping detail noise based on accumulated density causes some banding artefacts
		// so only use detail noise for the first two samples.
        bool use_detail_noise = i < 2;

        // Sample the cloud density at this point within the cone.
        density_along_cone += sample_cloud_density(p, height_fraction, float(i) * 0.5f, use_detail_noise);
	}

    // Get one more sample further away to account for shadows from distant clouds.
	_position += 32.0f * clouds.lightStepLength * _light_dir;

    // Compute height fraction for the distant position.
	float height_fraction = height_fraction_for_point(_position);

    // Sample the cloud density for the distant position.
	density_along_cone += sample_cloud_density(_position, height_fraction, 2.0f, false) * 3.0f;

	return density_along_cone;
}

// ------------------------------------------------------------------

float beer_lambert_law(float _density)
{
    return exp(-_density * clouds.precipitation);
}

// ------------------------------------------------------------------

float beer_law(float density)
{
	float d = -density * clouds.precipitation;
	return max(exp(d), exp(d * 0.5f) * 0.7f);
}

// ------------------------------------------------------------------

float henyey_greenstein_phase(float cos_angle, float g)
{
	float g2 = g * g;
	return ((1.0f - g2) / pow(1.0f + g2 - 2.0f * g * cos_angle, 1.5f)) / 4.0f * 3.1415f;
}

// ------------------------------------------------------------------

float powder_effect(float _density, float _cos_angle)
{
	float powder = 1.0f - exp(-_density * 2.0f);
	return mix(1.0f, powder, clamp((-_cos_angle * 0.5f) + 0.5f, 0.0f, 1.0f));
}

// ------------------------------------------------------------------

float calculate_light_energy(float _density, float _cos_angle, float _powder_density)
{
	float beer_powder = 2.0f * beer_law(_density) * powder_effect(_powder_density, _cos_angle);
	float HG = max(henyey_greenstein_phase(_cos_angle, clouds.henyeyGreensteinGForward), henyey_greenstein_phase(_cos_angle, clouds.henyeyGreensteinGBackward)) * 0.07f + 0.8f;
	return beer_powder * HG;
}

// ------------------------------------------------------------------

vec4 ray_march(vec3 _ray_origin, vec3 _ray_direction, float _cos_angle, float _step_size, float _num_steps)
{
	vec3  position            = _ray_origin;
	float step_increment      = 1.0f;
    float accum_transmittance = 1.0f;
    vec3  accum_scattering    = vec3(0.0f);
    float alpha               = 0.0f;

	vec3 sun_color = vec3(1.0);//GetSolarRadiance();

	for (float i = 0.0f; i < _num_steps; i+= step_increment)
	{
		float height_fraction    = height_fraction_for_point(position);
		float density            = sample_cloud_density(position, height_fraction, 0.0f, false);
        float step_transmittance = beer_lambert_law(density * _step_size);

        accum_transmittance *= step_transmittance;

		if (density > 0.0f)
		{
            alpha += (1.0f - step_transmittance) * (1.0f - alpha);

            float cone_density = sample_cloud_density_along_cone(position, vec3( 0.00, 0.84805, -0.52992));

            vec3 in_scattered_light = calculate_light_energy(cone_density * _step_size, _cos_angle, density * _step_size) * sun_color * clouds.sunLightFactor * alpha;
            vec3 ambient_light      = mix(clouds.cloudBaseColor, vec3(1.0), height_fraction) * clouds.ambientLightFactor; ///top colour

            accum_scattering += (ambient_light + in_scattered_light) * accum_transmittance * density;
		}

		position += _ray_direction * _step_size * step_increment;
	}

	return vec4(accum_scattering, alpha);
}

// ------------------------------------------------------------------
#endif

#ifdef ATMOSHPERE_VIEWER_IN_SPACE
RadianceSpectrum GetSkyRadiance(

    IN(TransmittanceTexture) transmittance_texture,
    IN(ReducedScatteringTexture) scattering_texture,
    IN(ReducedScatteringTexture) single_mie_scattering_texture,
    Position camera, IN(Direction) view_ray, Length shadow_length,
    IN(Direction) sun_direction, OUT(DimensionlessSpectrum) transmittance) {
  // Compute the distance to the top atmosphere boundary along the view ray,
  // assuming the viewer is in space (or NaN if the view ray does not intersect
  // the atmosphere).
  Length r = length(camera);
  Length rmu = dot(camera, view_ray);
  Length distance_to_top_atmosphere_boundary = -rmu -
      sqrt(rmu * rmu - r * r + top_radius * top_radius);
  // If the viewer is in space and the view ray intersects the  move
  // the viewer to the top atmosphere boundary (along the view ray):
  if (distance_to_top_atmosphere_boundary > 0.0 * m) {
    camera = camera + view_ray * distance_to_top_atmosphere_boundary;
    r = top_radius;
    rmu += distance_to_top_atmosphere_boundary;
  } else if (r > top_radius) {
    // If the view ray does not intersect the  simply return 0.
    transmittance = DimensionlessSpectrum(1.0);
    return RadianceSpectrum(0.0 * watt_per_square_meter_per_sr_per_nm);
  }
  // Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
  Number mu = rmu / r;
  Number mu_s = dot(camera, sun_direction) / r;
  Number nu = dot(view_ray, sun_direction);
  bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

  transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum(0.0) :
      GetTransmittanceToTopAtmosphereBoundary(
           transmittance_texture, r, mu);
  IrradianceSpectrum single_mie_scattering;
  IrradianceSpectrum scattering;
  if (shadow_length == 0.0 * m) {
    scattering = GetCombinedScattering(
         scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
  } else {
    // Case of light shafts (shadow_length is the total length noted l in our
    // paper): we omit the scattering between the camera and the point at
    // distance l, by implementing Eq. (18) of the paper (shadow_transmittance
    // is the T(x,x_s) term, scattering is the S|x_s=x+lv term).
    Length d = shadow_length;
    Length r_p =
        ClampRadius( sqrt(d * d + 2.0 * r * mu * d + r * r));
    Number mu_p = (r * mu + d) / r_p;
    Number mu_s_p = (r * mu_s + d * nu) / r_p;

    scattering = GetCombinedScattering(
         scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
    DimensionlessSpectrum shadow_transmittance =
        GetTransmittance( transmittance_texture,
            r, mu, shadow_length, ray_r_mu_intersects_ground);
    scattering = scattering * shadow_transmittance;
    single_mie_scattering = single_mie_scattering * shadow_transmittance;
  }
  return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
      MiePhaseFunction(mie_phase_function_g, nu);
}
#else
RadianceSpectrum GetSkyRadiance(

    IN(TransmittanceTexture) transmittance_texture,
    IN(ReducedScatteringTexture) scattering_texture,
    IN(ReducedScatteringTexture) single_mie_scattering_texture,
    IN(Length) rmu, IN(Length) nu, IN(Direction) view_ray, Length shadow_length,
    IN(Direction) sun_direction, OUT(DimensionlessSpectrum) transmittance) {
  // Compute the distance to the top atmosphere boundary along the view ray,
  // assuming the viewer is in space (or NaN if the view ray does not intersect
  // the atmosphere).
  Length r = positional.cameraPosR.w;
  Length distance_to_top_atmosphere_boundary = -rmu -
      sqrt(rmu * rmu - r * r + top_radius * top_radius);

  // Compute the r, mu, mu_s and nu parameters needed for the texture lookups.
  Number mu = rmu / r;
  Number mu_s = positional.mu_s;
  bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

  transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum(0.0) :
      GetTransmittanceToTopAtmosphereBoundary(
           transmittance_texture, r, mu);
  IrradianceSpectrum single_mie_scattering;
  IrradianceSpectrum scattering;
  if (shadow_length == 0.0 * m) {
    scattering = GetCombinedScattering(
         scattering_texture, single_mie_scattering_texture,
        r, mu, mu_s, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
  } else {
    // Case of light shafts (shadow_length is the total length noted l in our
    // paper): we omit the scattering between the camera and the point at
    // distance l, by implementing Eq. (18) of the paper (shadow_transmittance
    // is the T(x,x_s) term, scattering is the S|x_s=x+lv term).
    Length d = shadow_length;
    Length r_p =
        ClampRadius( sqrt(d * d + 2.0 * r * mu * d + r * r));
    Number mu_p = (r * mu + d) / r_p;
    Number mu_s_p = (r * mu_s + d * nu) / r_p;

    scattering = GetCombinedScattering(
         scattering_texture, single_mie_scattering_texture,
        r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
        single_mie_scattering);
    DimensionlessSpectrum shadow_transmittance =
        GetTransmittance( transmittance_texture,
            r, mu, shadow_length, ray_r_mu_intersects_ground);
    scattering = scattering * shadow_transmittance;
    single_mie_scattering = single_mie_scattering * shadow_transmittance;
  }
  return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
      MiePhaseFunction(mie_phase_function_g, nu);
}
#endif

void main()
{
    // Generate a camera ray to the current fragment.
	vec3 view_direction = normalize(inRay);

    float r = positional.cameraPosR.w;
    float rmu = dot(positional.cameraPosR.xyz, view_direction);
    Number nu = dot(view_direction, positional.sunDirectionExp.xyz);

    float cloudsMin = bottom_radius + clouds.cloudMinHeight;
    float cloudsMax = bottom_radius + clouds.cloudMaxHeight;

    float distance_to_top_clouds_boundary = rmu +
      sqrt(rmu * rmu - r * r + top_radius * top_radius);
    float distance_to_bottom_clouds_boundary = rmu +
      sqrt(rmu * rmu - r * r + bottom_radius * bottom_radius);

    // Compute the radiance of the sky.
	vec3 transmittance;
	vec3 radiance = GetSkyRadiance(s_Transmittance,
                                   s_Scattering,
                                s_SingleMieScattering, rmu, nu,
                                view_direction,
                                4.0,
                                positional.sunDirectionExp.xyz,
                                transmittance) * SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;

	// If the view ray intersects the Sun, add the Sun radiance.
	if (dot(view_direction, positional.sunDirectionExp.xyz) > settings.sunSize.y)
		radiance = radiance + transmittance * GetSolarRadiance();

	vec4 sky = vec4(pow(vec3(1.0) - exp(-radiance / settings.whitePoint.rbg * positional.sunDirectionExp.a), vec3(1.0 / 2.2)), 1.0);
/*
    vec4 cloud = vec4(0.0);
    if(distance_to_bottom_clouds_boundary < distance_to_top_clouds_boundary)
    {
        // Get a random number that we'll use to jitter our ray.
        const float rng = blue_noise();

        // The maximum number of ray march steps to use.
        const float max_steps = clouds.maxNumSteps;

        // The minimum number of ray march steps to use with an added offset to prevent banding.
        const float min_steps = (clouds.maxNumSteps * 0.5f) + (rng * 2.0f);

        // The number of ray march steps is determined depending on how steep the viewing angle is.
        float num_steps = mix(max_steps, min_steps, view_direction.y);

        // Using the number of steps we can determine the step size of our ray march.
        float step_size = (distance_to_top_clouds_boundary - distance_to_bottom_clouds_boundary) / num_steps;

        // Jitter the ray to prevent banding.
        vec3 ray_start = positional.cameraPosR.xyz + (distance_to_bottom_clouds_boundary * view_direction);// + step_size * view_direction * rng;

        cloud = ray_march(ray_start, view_direction, nu, step_size, num_steps);
    }
    //vec3 cameraPos = pc.modelView[3].xyz;
*/

// Generate a camera ray to the current fragment.
	Ray ray = generate_ray();
    //ray.direction = view_direction;

	// Figure out where our ray will intersect the sphere that represents the beginning of the cloud layer.
	vec3 ray_start = ray_sphere_intersection(ray, vec3(0.0), bottom_radius);

	// Similarly, figure out where the ray intersects the end of the cloud layer.
	vec3 ray_end   = ray_sphere_intersection(ray, vec3(0.0), top_radius);

	// Get a random number that we'll use to jitter our ray.
	const float rng = blue_noise();

	// The maximum number of ray march steps to use.
	const float max_steps = clouds.maxNumSteps;

	// The minimum number of ray march steps to use with an added offset to prevent banding.
	const float min_steps = (clouds.maxNumSteps * 0.5f) + (rng * 2.0f);

	// The number of ray march steps is determined depending on how steep the viewing angle is.
	float num_steps = mix(max_steps, min_steps, ray.direction.y);

	// Using the number of steps we can determine the step size of our ray march.
	float step_size = length(ray_start - ray_end) / num_steps;

	// Jitter the ray to prevent banding.
	ray_start += step_size * ray.direction * rng;

	float cos_angle = dot(ray.direction, vec3( 0.00, 0.84805, -0.52992));
	vec4  cloud    = ray_march(ray_start, ray.direction, cos_angle, step_size, num_steps);

    //outColor = vec4(clouds.rgb + (1.0f - clouds.a) * sky.rgb, 1.0);
    outColor = vec4(cloud.rgb, 1.0);
}
