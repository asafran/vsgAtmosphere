#version 450

#include "compute_uniforms.glsl"
#include "constants.glsl"
#include "functions.glsl"

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout (local_size_x_id = 0, local_size_y_id = 0, local_size_z = 1) in;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

layout (set = 0, binding = 0, rgba32f) uniform image2D transmittance;

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    ivec2 coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 frag_coord = coord + vec2(0.5, 0.5);
    imageStore(transmittance, coord, vec4(ComputeTransmittanceToTopAtmosphereBoundaryTexture(frag_coord), 1.0));
}

// ------------------------------------------------------------------
