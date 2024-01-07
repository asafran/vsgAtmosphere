#version 450
#extension GL_ARB_separate_shader_objects : enable
#pragma import_defines (ATMOSHPERE_RADIANCE, ATMOSHPERE_VIEWER_IN_SPACE)

#include "rendering_constants.glsl"
#include "constants.glsl"
#include "functions.glsl"


layout(location = 0) in vec3 inRay;
layout(location = 1) in vec2 inUV;
layout(location = 0) out vec4 outColor;

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

layout(push_constant) uniform PushConstants {
    mat4 projection;
    mat4 modelView;
} pc; //inverse

// ------------------------------------------------------------------
// FUNCTIONS --------------------------------------------------------
// ------------------------------------------------------------------


Luminance3 GetSolarRadiance()
{
  return solar_irradiance / (PI * sun_angular_radius * sun_angular_radius);
}

#ifdef ATMOSHPERE_VIEWER_IN_SPACE
RadianceSpectrum GetSkyRadiance(IN(Direction) view_ray, IN(Number) nu, OUT(DimensionlessSpectrum) transmittance)
{
  // Compute the distance to the top atmosphere boundary along the view ray,
  // assuming the viewer is in space (or NaN if the view ray does not intersect
  // the atmosphere).
  Length r = positional.radius;
  Position camera = positional.cameraPos;
  Number rmu = dot(positional.cameraPos, view_ray);

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
  Number mu_s = dot(camera, positional.sunDirection) / r;
  bool ray_r_mu_intersects_ground = RayIntersectsGround( r, mu);

  transmittance = ray_r_mu_intersects_ground ? DimensionlessSpectrum(0.0) :
      GetTransmittanceToTopAtmosphereBoundary(
           transmittance_texture, r, mu);
  IrradianceSpectrum single_mie_scattering;
   // Case of light shafts (shadow_length is the total length noted l in our
    // paper): we omit the scattering between the camera and the point at
    // distance l, by implementing Eq. (18) of the paper (shadow_transmittance
    // is the T(x,x_s) term, scattering is the S|x_s=x+lv term).
  Length d = 4.0;
  Length r_p =
      ClampRadius(sqrt(d * d + 2.0 * r * mu * d + r * r));
  Number mu_p = (r * mu + d) / r_p;
  Number mu_s_p = (r * mu_s + d * nu) / r_p;

  IrradianceSpectrum scattering = GetCombinedScattering(
    scattering_texture, single_mie_scattering_texture,
      r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground,
      single_mie_scattering);
  DimensionlessSpectrum shadow_transmittance =
        GetTransmittance(transmittance_texture,
            r, mu, d, ray_r_mu_intersects_ground);
    scattering = scattering * shadow_transmittance;
    single_mie_scattering = single_mie_scattering * shadow_transmittance;

  return scattering * RayleighPhaseFunction(nu) + single_mie_scattering *
      MiePhaseFunction(mie_phase_function_g, nu);
}
#else
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
#endif

void main()
{
	vec3 view_direction = normalize(inRay);
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

	outColor = vec4(pow(vec3(1.0) - exp(-radiance / settings.whitePoint.rbg * positional.exposure), vec3(1.0 / 2.2)), 1.0);
}
