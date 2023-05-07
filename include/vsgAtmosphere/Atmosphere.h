#ifndef ATMOSPHEREMODEL_H
#define ATMOSPHEREMODEL_H

#include <string>
#include <vector>
#include <vsg/app/EllipsoidModel.h>
#include <vsg/app/Window.h>
#include <vsg/app/WindowResizeHandler.h>
#include <vsg/core/Data.h>
#include <vsg/state/ImageInfo.h>
#include <vsg/state/ShaderStage.h>
#include <vsg/utils/ShaderSet.h>

namespace atmosphere {

struct DensityProfileLayer
{
    DensityProfileLayer(double in_width, double in_exp_term, double in_exp_scale, double in_linear_term, double in_constant_term, double lengthUnitInMeters)
    {
        width = static_cast<float>(in_width / lengthUnitInMeters);
        exp_term = static_cast<float>(in_exp_term);
        exp_scale = static_cast<float>(in_exp_scale * lengthUnitInMeters);
        linear_term = static_cast<float>(in_linear_term * lengthUnitInMeters);
        constant_term = static_cast<float>(in_constant_term);
    }

    DensityProfileLayer() {}

    float width = 0.0f;
    float exp_term = 0.0f;
    float exp_scale = 0.0f;
    float linear_term = 0.0f;
    float constant_term = 0.0f;
};

struct RuntimeSettings
{
    vsg::vec4 whitePointExp;
    vsg::vec2 sunSize;
};

class AtmosphereModel : public vsg::Inherit<vsg::Object, AtmosphereModel>
{
public:

    int numThreads = 8;
    int numViewerThreads = 32;

    int transmittanceWidth = 256;
    int transmittanceHeight = 64;

    int scaterringR = 32;
    int scaterringMU = 128;
    int scaterringMU_S = 32;
    int scaterringNU = 8;

    int irradianceWidth = 64;
    int irradianceHeight = 16;

    int cubeSize = 1024;

    /// <summary>
    /// The wavelength values, in nanometers, and sorted in increasing order, for
    /// which the solar_irradiance, rayleigh_scattering, mie_scattering,
    /// mie_extinction and ground_albedo samples are provided. If your shaders
    /// use luminance values (as opposed to radiance values, see above), use a
    /// large number of wavelengths (e.g. between 15 and 50) to get accurate
    /// results (this number of wavelengths has absolutely no impact on the
    /// shader performance).
    /// </summary>
    std::vector<double> waveLengths;

    /// <summary>
    /// The solar irradiance at the top of the atmosphere, in W/m^2/nm. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> solarIrradiance;

    /// <summary>
    /// The sun's angular radius, in radians. Warning: the implementation uses
    /// approximations that are valid only if this value is smaller than 0.1.
    /// </summary>
    float sunAngularRadius;

    /// <summary>
    /// The distance between the planet center and the bottom of the atmosphere in m.
    /// </summary>
    double bottomRadius;

    /// <summary>
    /// The distance between the planet center and the top of the atmosphere in m.
    /// </summary>
    double topRadius;

    /// <summary>
    /// The density profile of air molecules, i.e. a function from altitude to
    /// dimensionless values between 0 (null density) and 1 (maximum density).
    /// Layers must be sorted from bottom to top. The width of the last layer is
    /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
    /// layers can be specified.
    /// </summary>
    DensityProfileLayer rayleighDensityLayer;

    /// <summary>
    /// The scattering coefficient of air molecules at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
    /// to 'rayleigh_scattering' times 'rayleigh_density' at this altitude. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> rayleighScattering;

    /// <summary>
    /// The density profile of aerosols, i.e. a function from altitude to
    /// dimensionless values between 0 (null density) and 1 (maximum density).
    /// Layers must be sorted from bottom to top. The width of the last layer is
    /// ignored, i.e. it always extend to the top atmosphere boundary. At most 2
    /// layers can be specified.
    /// </summary>
    DensityProfileLayer mieDensityLayer;

    /// <summary>
    /// The scattering coefficient of aerosols at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The scattering coefficient at altitude h is equal
    /// to 'mie_scattering' times 'mie_density' at this altitude. This vector
    /// must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> mieScattering;

    /// <summary>
    /// The extinction coefficient of aerosols at the altitude where their
    /// density is maximum (usually the bottom of the atmosphere), as a function
    /// of wavelength, in m^-1. The extinction coefficient at altitude h is equal
    /// to 'mie_extinction' times 'mie_density' at this altitude. This vector
    /// must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> mieExtinction;

    /// <summary>
    /// The asymetry parameter for the Cornette-Shanks phase function for the aerosols.
    /// </summary>
    float miePhaseFunction_g;

    /// <summary>
    /// The density profile of air molecules that absorb light (e.g. ozone), i.e.
    /// a function from altitude to dimensionless values between 0 (null density)
    /// and 1 (maximum density). Layers must be sorted from bottom to top. The
    /// width of the last layer is ignored, i.e. it always extend to the top
    /// atmosphere boundary. At most 2 layers can be specified.
    /// </summary>
    DensityProfileLayer absorptionDensityLayer0;
    DensityProfileLayer absorptionDensityLayer1;

    /// <summary>
    /// The extinction coefficient of molecules that absorb light (e.g. ozone) at
    /// the altitude where their density is maximum, as a function of wavelength,
    /// in m^-1. The extinction coefficient at altitude h is equal to
    /// 'absorption_extinction' times 'absorption_density' at this altitude. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> absorptionExtinction;

    /// <summary>
    /// The average albedo of the ground, as a function of wavelength. This
    /// vector must have the same size as the wavelengths parameter.
    /// </summary>
    std::vector<double> groundAlbedo;

    /// <summary>
    /// The maximum Sun zenith angle for which atmospheric scattering must be
    /// precomputed, in radians (for maximum precision, use the smallest Sun
    /// zenith angle yielding negligible sky light radiance values. For instance,
    /// for the Earth case, 102 degrees is a good choice for most cases (120
    /// degrees is necessary for very high exposure values).
    /// </summary>
    double maxSunZenithAngle;

    /// <summary>
    /// The length unit used in your shaders and meshes. This is the length unit
    /// which must be used when calling the atmosphere model shader functions.
    /// </summary>
    double lengthUnitInMeters;

    vsg::ref_ptr<vsg::Data> transmittanceData;
    vsg::ref_ptr<vsg::Data> scatteringData;
    vsg::ref_ptr<vsg::Data> irradianceData;
    vsg::ref_ptr<vsg::Data> mieScatteringData;

    vsg::ref_ptr<vsg::DescriptorSetLayout> viewDescriptorSetLayout;

    vsg::ref_ptr<vsg::ShaderCompileSettings> compileSettings;

    vsg::ref_ptr<vsg::Value<atmosphere::RuntimeSettings>> runtimeSettings;

    vsg::dvec3 sunDirection = {0.0, std::sin(vsg::PI), std::cos(vsg::PI)};

    int precomputedWavelenghts = 15;

private:

    struct Parameters
    {
        vsg::vec4 solar_irradiance;
        vsg::vec4 rayleigh_scattering;
        vsg::vec4 mie_scattering;
        vsg::vec4 mie_extinction;
        vsg::vec4 ground_albedo;
        vsg::vec4 absorption_extinction;
    };

    vsg::ref_ptr<vsg::Device> _device;
    vsg::ref_ptr<vsg::PhysicalDevice> _physicalDevice;
    vsg::ref_ptr<vsg::Options> _options;

    vsg::ref_ptr<vsg::ImageInfo> _transmittanceTexture;

    vsg::ref_ptr<vsg::ImageInfo> _irradianceTexture;
    vsg::ref_ptr<vsg::ImageInfo> _deltaIrradianceTexture;

    vsg::ref_ptr<vsg::ImageInfo> _deltaRayleighScatteringTexture;
    vsg::ref_ptr<vsg::ImageInfo> _deltaMieScatteringTexture;
    vsg::ref_ptr<vsg::ImageInfo> _scatteringTexture;
    vsg::ref_ptr<vsg::ImageInfo> _singleMieScatteringTexture;

    vsg::ref_ptr<vsg::ImageInfo> _deltaScatteringDensityTexture;
    vsg::ref_ptr<vsg::ImageInfo> _deltaMultipleScatteringTexture;

    vsg::ref_ptr<vsg::ImageInfo> _cubeMap;

    vsg::ref_ptr<vsg::BufferInfo> _directLight;
    vsg::ref_ptr<vsg::BufferInfo> _ambientLight;

    vsg::ShaderStage::SpecializationConstants _computeConstants;
    vsg::ShaderStage::SpecializationConstants _renderConstants;

public:
    AtmosphereModel(vsg::ref_ptr<vsg::Device> device, vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice, vsg::ref_ptr<vsg::Options> options);
    ~AtmosphereModel();

    void initialize(int scatteringOrders);

    vsg::vec3 convertSpectrumToLinearSrgb(double c);

    vsg::ref_ptr<vsg::CommandGraph> createCubeMapGraph(vsg::ref_ptr<vsg::vec4Value> camera);
    vsg::ref_ptr<vsg::Node> createEnvMap();
    vsg::ref_ptr<vsg::Node> createSky();
    vsg::ref_ptr<vsg::Node> createSkyView(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Camera> camera);
    vsg::ref_ptr<vsg::Node> createSkyBox();
    vsg::ref_ptr<vsg::ShaderSet> phongShaderSet();

    void setSunAngle(double radians);
    void setExposure(double value);

private:
    void generateTextures();

    vsg::ref_ptr<vsg::ImageInfo> generate2D(uint32_t width, uint32_t height, bool init = false);
    vsg::ref_ptr<vsg::ImageInfo> generate3D(uint32_t width, uint32_t height, uint32_t depth, bool init = false);
    vsg::ref_ptr<vsg::ImageInfo> generateCubemap(uint32_t size);
    //vsg::ref_ptr<vsg::BufferInfo> generateVec4Array(uint32_t size);

    int scatteringWidth() const { return scaterringNU * scaterringMU_S; }
    int scatteringHeight() const { return scaterringMU; }
    int scatteringDepth() const { return scaterringR; }

    vsg::ref_ptr<vsg::BindComputePipeline> bindCompute(const vsg::Path& filename, vsg::ref_ptr<vsg::PipelineLayout> pipelineLayout) const;
    vsg::ref_ptr<vsg::DescriptorSet> bindTransmittance() const;
    vsg::ref_ptr<vsg::DescriptorSet> bindDirectIrradiance() const;
    vsg::ref_ptr<vsg::DescriptorSet> bindSingleScattering() const;

    vsg::ref_ptr<vsg::DescriptorSet> bindScatteringDensity() const;
    vsg::ref_ptr<vsg::DescriptorSet> bindIndirectIrradiance() const;
    vsg::ref_ptr<vsg::DescriptorSet> bindMultipleScattering() const;

    void computeParameters(vsg::ref_ptr<vsg::Value<Parameters>> parameters, const vsg::vec3 &lambdas) const;

    vsg::ref_ptr<vsg::DescriptorSetLayout> parametersLayout() const;
    vsg::ref_ptr<vsg::DescriptorSetLayout> orderLayout() const;

    void assignComputeConstants();
    void assignRenderConstants();

    vsg::ref_ptr<vsg::Data> mapData(vsg::ref_ptr<vsg::ImageView> view, uint32_t width, uint32_t height);
    vsg::ref_ptr<vsg::Data> mapData(vsg::ref_ptr<vsg::ImageView> view, uint32_t width, uint32_t height, uint32_t depth);

    float exposure = 10.0f;

    friend class AtmosphereLighting;
};

extern vsg::ref_ptr<AtmosphereModel> createAtmosphereModel(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::EllipsoidModel> eps, vsg::ref_ptr<vsg::Options> options);

}
#endif // ATMOSPHEREMODEL_H
