#version 450

#include "compute_uniforms.glsl"
#include "constants.glsl"
#include "functions.glsl"

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x_id = 0, local_size_y_id = 0, local_size_z = 1) in;

// ------------------------------------------------------------------
// IMAGES -----------------------------------------------------------
// ------------------------------------------------------------------

layout (set = 0, binding = 0, rgba32f) uniform image2D delta_irradiance;
layout (set = 0, binding = 1, rgba32f) uniform image2D irradiance;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

//uniform vec4 blend;
layout (set = 2, binding = 0) uniform ScatteringOrder
{
    int scattering_order;
};

layout (set = 0, binding = 2) uniform sampler3D single_rayleigh_scattering;
layout (set = 0, binding = 3) uniform sampler3D single_mie_scattering;
layout (set = 0, binding = 4) uniform sampler3D multiple_scattering;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec2 coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 frag_coord = coord + vec2(0.5, 0.5);

    vec3 delta_irradiance_value = ComputeIndirectIrradianceTexture(single_rayleigh_scattering, single_mie_scattering, multiple_scattering, frag_coord, scattering_order - 1);

    imageStore(delta_irradiance, coord, vec4(delta_irradiance_value, 1.0));

    imageStore(irradiance, coord, imageLoad(irradiance, coord) + (luminance_from_radiance * vec4(delta_irradiance_value, 0.0)));
}

// ------------------------------------------------------------------
