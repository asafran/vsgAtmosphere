layout (set = 1, binding = 0) uniform LuminanceFromRadiance
{
    mat4 luminance_from_radiance;
};

layout (set = 1, binding = 1) uniform Parameters
{
    vec4 solar_irradiance;
    vec4 rayleigh_scattering;
    vec4 mie_scattering;
    vec4 mie_extinction;
    vec4 ground_albedo;
    vec4 absorption_extinction;
};

struct DensityProfileLayer
{
  float width;
  float exp_term;
  float exp_scale;
  float linear_term;
  float constant_term;
};

layout (std430, set = 1, binding = 2) buffer DensityProfileLayers
{
    DensityProfileLayer[4] layers;
};
