layout (constant_id = 18) const float rsX = 0.0f;
layout (constant_id = 19) const float rsY = 0.0f;
layout (constant_id = 20) const float rsZ = 0.0f;
const vec3 rayleigh_scattering = vec3(rsX, rsY, rsZ);

layout (constant_id = 21) const float mieX = 0.0f;
layout (constant_id = 22) const float mieY = 0.0f;
layout (constant_id = 23) const float mieZ = 0.0f;
const vec3 mie_scattering = vec3(mieX, mieY, mieZ);

layout (constant_id = 27) const float irrX = 0.0f;
layout (constant_id = 28) const float irrY = 0.0f;
layout (constant_id = 29) const float irrZ = 0.0f;
const vec3 solar_irradiance = vec3(irrX, irrY, irrZ);


layout (constant_id = 24) const float aX = 0.0f;
layout (constant_id = 25) const float aY = 0.0f;
layout (constant_id = 26) const float aZ = 0.0f;
const vec3 absorption_extinction = vec3(aX, aY, aZ);

layout (constant_id = 30) const float mieexX = 0.0f;
layout (constant_id = 31) const float mieexY = 0.0f;
layout (constant_id = 32) const float mieexZ = 0.0f;
const vec3 mie_extinction = vec3(mieexX, mieexY, mieexZ);

layout (constant_id = 33) const float gX = 0.0f;
layout (constant_id = 34) const float gY = 0.0f;
layout (constant_id = 35) const float gZ = 0.0f;
const vec3 ground_albedo = vec3(gX, gY, gZ);


layout (constant_id = 41) const float skyX = 0.0f;
layout (constant_id = 42) const float skyY = 0.0f;
layout (constant_id = 43) const float skyZ = 0.0f;
const vec3 SKY_SPECTRAL_RADIANCE_TO_LUMINANCE = vec3(skyX, skyY, skyZ);

layout (constant_id = 44) const float sunX = 0.0f;
layout (constant_id = 45) const float sunY = 0.0f;
layout (constant_id = 46) const float sunZ = 0.0f;
const vec3 SUN_SPECTRAL_RADIANCE_TO_LUMINANCE = vec3(sunX, sunY, sunZ);

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
