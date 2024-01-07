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

layout (set = 0, binding = 0, rgba32f) uniform image3D delta_scattering_density;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

//uniform vec4 blend;
//uniform int layer;

layout (set = 0, binding = 1) uniform sampler2D transmittance;
layout (set = 0, binding = 2) uniform sampler3D single_rayleigh_scattering;
layout (set = 0, binding = 3) uniform sampler3D single_mie_scattering;
layout (set = 0, binding = 4) uniform sampler3D multiple_scattering;
layout (set = 0, binding = 5) uniform sampler2D irradiance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    //coord.z = layer;
    vec3 frag_coord = coord + vec3(0.5,0.5,0.5);

    vec3 scattering_density = ComputeScatteringDensityTexture(transmittance, single_rayleigh_scattering, single_mie_scattering, multiple_scattering, irradiance, frag_coord, scattering_order);

    imageStore(delta_scattering_density, coord, vec4(scattering_density, 1.0));
}

// ------------------------------------------------------------------
