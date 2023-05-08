#version 450
#extension GL_ARB_separate_shader_objects : enable

#include "rendering_constants.glsl"
#include "constants.glsl"
#include "functions.glsl"

// ------------------------------------------------------------------
// CONSTANTS --------------------------------------------------------
// ------------------------------------------------------------------
layout (local_size_x_id = 1, local_size_y_id = 1, local_size_z = 1) in;
layout (constant_id = 2) const int CUBE_SIZE = 1024;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

layout(set = 0, binding = 0, std140) uniform Settings
{
	vec4 white_point_exp;
	vec4 sun_direction;
	vec2 sun_size;
} settings;

layout(set = 0, binding = 1) uniform sampler2D transmittance_texture;
layout(set = 0, binding = 2) uniform sampler2D irradiance_texture;
layout(set = 0, binding = 3) uniform sampler3D scattering_texture;
layout(set = 0, binding = 4) uniform sampler3D single_mie_scattering_texture;

layout(set = 0, binding = 5, rgba32f) uniform imageCube cubemap;

layout(push_constant) uniform PushConstants {
    vec4 position;
} camera;

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
      return GetSkyRadiance(transmittance_texture,
          scattering_texture, single_mie_scattering_texture,
          camera, view_ray, shadow_length, sun_direction, transmittance) *
          SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
    Luminance3 GetSkyRadianceToPoint(
        Position camera, Position point, Length shadow_length,
        Direction sun_direction, out DimensionlessSpectrum transmittance) {
      return GetSkyRadianceToPoint(transmittance_texture,
          scattering_texture, single_mie_scattering_texture,
          camera, point, shadow_length, sun_direction, transmittance) *
          SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }
    Illuminance3 GetSunAndSkyIrradiance(
       Position p, Direction normal, Direction sun_direction,
       out IrradianceSpectrum sky_irradiance) {
      IrradianceSpectrum sun_irradiance = GetSunAndSkyIrradiance(
          transmittance_texture, irradiance_texture, p, normal,
          sun_direction, sky_irradiance);
      sky_irradiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
      return sun_irradiance * SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
    }

/*
RadianceSpectrum GetSolarRadiance() {
		return solar_irradiance / (PI * sun_angular_radius * sun_angular_radius);
	}

// ------------------------------------------------------------------

RadianceSpectrum GetSkyRadiance(
		Position camera, Direction view_ray, Length shadow_length,
		Direction sun_direction, out DimensionlessSpectrum transmittance) {
		return GetSkyRadiance(transmittance_texture,
			scattering_texture, single_mie_scattering_texture,
			camera, view_ray, shadow_length, sun_direction, transmittance);
	}
*/
// ------------------------------------------------------------------

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

	// Compute the radiance of the sky.
	vec3 transmittance;
	vec3 radiance = GetSkyRadiance(camera.position.xyz, view_direction, 4.0, settings.sun_direction.xyz, transmittance);

	// If the view ray intersects the Sun, add the Sun radiance.
	if (dot(view_direction, settings.sun_direction.xyz) > settings.sun_size.y)
		radiance = radiance + transmittance * GetSolarRadiance();

	ivec3 texel = ivec3(gl_GlobalInvocationID);
		
	imageStore(cubemap, texel, vec4(pow(vec3(1.0) - exp(-radiance / settings.white_point_exp.rbg * settings.white_point_exp.a), vec3(1.0 / 2.2)), 1.0));
}
