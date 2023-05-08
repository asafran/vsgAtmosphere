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

layout (set = 0, binding = 0, rgba32f) uniform image3D delta_rayleigh_scattering;
layout (set = 0, binding = 1, rgba32f) uniform image3D delta_mie_scattering;
layout (set = 0, binding = 2, rgba32f) uniform image3D scattering;
layout (set = 0, binding = 3, rgba32f) uniform image3D single_mie_scattering;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

//uniform vec4 blend;
//uniform int layer;
layout (set = 0, binding = 4) uniform sampler2D transmittance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec3 coord = ivec3(gl_GlobalInvocationID);
    //coord.z = layer;
    vec3 frag_coord = coord + vec3(0.5,0.5,0.5);

    vec3 delta_rayleigh, delta_mie;
    ComputeSingleScatteringTexture(transmittance, frag_coord, delta_rayleigh, delta_mie);

    vec4 delta_rayleigh_rgba = vec4(delta_rayleigh, 1);
    vec4 delta_mie_rgba = vec4(delta_mie, 1);

    imageStore(delta_rayleigh_scattering, coord, delta_rayleigh_rgba);

    imageStore(delta_mie_scattering, coord, delta_mie_rgba);

    vec4 scattering_source = imageLoad(scattering, coord);

    vec4 single_mie_scattering_source = imageLoad(single_mie_scattering, coord);

    imageStore(scattering, coord, scattering_source + vec4((luminance_from_radiance * delta_rayleigh_rgba).rgb, (luminance_from_radiance * delta_mie_rgba).r));

    imageStore(single_mie_scattering, coord, single_mie_scattering_source + (luminance_from_radiance * delta_mie_rgba));
}

// ------------------------------------------------------------------
