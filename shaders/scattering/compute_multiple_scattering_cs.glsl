#version 450

#include "compute_uniforms.glsl"
#include "constants.glsl"
#include "functions.glsl"

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x_id = 0, local_size_y_id = 0, local_size_z_id = 0) in;

// ------------------------------------------------------------------
// IMAGES -----------------------------------------------------------
// ------------------------------------------------------------------

layout (set = 0, binding = 0, rgba32f) uniform image3D delta_multiple_scattering;
layout (set = 0, binding = 1, rgba32f) uniform image3D scattering;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

//uniform vec4 blend;
//uniform int layer;
layout (set = 2, binding = 0) uniform ScatteringOrder
{
    int scattering_order;
};

layout (set = 0, binding = 2) uniform sampler2D transmittance;
layout (set = 0, binding = 3) uniform sampler3D delta_scattering_density;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    //coord.z = layer;

    vec3 frag_coord = coord + vec3(0.5, 0.5, 0.5);

    float nu;
    vec4 delta_multiple_scattering_value = vec4(ComputeMultipleScatteringTexture(transmittance, delta_scattering_density, frag_coord, nu), 1.0);

    imageStore(delta_multiple_scattering, coord, delta_multiple_scattering_value);

    imageStore(scattering, coord, imageLoad(scattering, coord) + vec4((luminance_from_radiance * delta_multiple_scattering_value).rgb / RayleighPhaseFunction(nu), 0.0));
}

// ------------------------------------------------------------------
