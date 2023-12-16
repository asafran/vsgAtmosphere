#include "noise.glsl"

// ------------------------------------------------------------------
// INPUTS -----------------------------------------------------------
// ------------------------------------------------------------------

layout(local_size_x_id = 0, local_size_y_id = 0, local_size_z_id = 0) in;

// ------------------------------------------------------------------
// UNIFORMS ---------------------------------------------------------
// ------------------------------------------------------------------

<<<<<<< Updated upstream
layout(binding = 0, rgba16f) uniform image3D detail_noise;

layout (constant_id = 41) const int size = 32;
=======
layout(binding = 1, rgba16f) uniform image3D detail_noise;

layout (constant_id = 1) const int size = 32;
>>>>>>> Stashed changes

// ------------------------------------------------------------------
// MAIN -------------------------------------------------------------
// ------------------------------------------------------------------

void main()
{
    vec3 tex_coord = (vec3(gl_GlobalInvocationID) + vec3(0.5f)) / float(size);

    float freq = 8.0f;

    float worley0 = worley_fbm(tex_coord, freq);
    float worley1 = worley_fbm(tex_coord, freq * 2.0f);
    float worley2 = worley_fbm(tex_coord, freq * 4.0f);
 
    vec4 worley = vec4(worley0, worley1, worley2, 0.0f); 
    
    imageStore(detail_noise, ivec3(gl_GlobalInvocationID), worley);
}

// ------------------------------------------------------------------
