
layout(location = 0) in vec3 inRay;
layout(location = 1) in vec2 inUV;
layout(location = 0) out vec4 outColor;

layout(push_constant) uniform PushConstants {
    mat4 projection;
    mat4 modelView;
} pc; //inverse

layout(set = 1, binding = 2) uniform sampler2D transmittance_texture;
layout(set = 1, binding = 3) uniform sampler2D irradiance_texture;
layout(set = 1, binding = 4) uniform sampler3D scattering_texture;
layout(set = 1, binding = 5) uniform sampler3D single_mie_scattering_texture;

const float time = 0.0;
const float cirrus = 0.4;
const float cumulus = 0.8;

layout(set = 1, binding = 6, std140) uniform Settings
{
	vec4 whitePoint;
	vec2 sunSize;
} settings;

layout(set = 1, binding = 7, std140) uniform Positional
{
	vec4 sunDirectionExp;
    vec4 cameraPos;
} positional;

layout(set = 1, binding = 8) uniform sampler2D weather_texture;
layout(set = 1, binding = 9) uniform sampler3D shapenoise_texture;
layout(set = 1, binding = 10) uniform sampler3D detailnoise_texture;
layout(set = 1, binding = 11) uniform sampler2D bluenoise_texture;

layout(set = 1, binding = 12, std140) uniform CloudProperties{
    vec3 density_to_sigma_s;
    float phase_g;
    vec3 density_to_sigma_t;
    int primary_ray_marching_steps;
    vec3 box_min;
    int secondary_ray_marching_steps;
    vec3 box_max;
    int enable_multi_scattering;
    float g_c;
    float g_d;
    float wc0;
    float wc1;
    float wh;
    float shape_tile;
    float detail_tile;
    float blend_alpha;
};

layout(set = 1, binding = 13, std140) uniform CloudParams{
    vec3 camera_dir;
    float scale;//tan(fov/2)
    vec3 up;
    float w_over_h;
    vec3 right;
    vec3 view_pos;
    vec3 sun_dir;
};
