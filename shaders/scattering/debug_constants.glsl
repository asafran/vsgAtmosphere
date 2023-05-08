const vec3 rayleigh_scattering = vec3(0.0058, 0.01356, 0.0331);

const vec3 mie_scattering = vec3(0.004, 0.004, 0.004);

const vec3 solar_irradiance = vec3(1.474, 1.8504, 1.91198);

const vec3 absorption_extinction = vec3(0.00065, 0.00188, 0.00009);

const vec3 mie_extinction = vec3(0.00444, 0.00444, 0.00444);

const vec3 ground_albedo = vec3(0.1, 0.1, 0.1);

const vec3 SKY_SPECTRAL_RADIANCE_TO_LUMINANCE = vec3(114974.91406, 71305.95313, 65310.54688);

const vec3 SUN_SPECTRAL_RADIANCE_TO_LUMINANCE = vec3(98242.78906, 69954.39844, 66475.01563);

struct DensityProfileLayer
{
  float width;
  float exp_term;
  float exp_scale;
  float linear_term;
  float constant_term;
};

const DensityProfileLayer[4] layers = DensityProfileLayer[4](DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0),
                                                             DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0),
                                                             DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0),
                                                             DensityProfileLayer(0.0, 0.0, 0.0, 0.0, 0.0));
