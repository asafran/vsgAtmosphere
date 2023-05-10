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

#include <vsg/io/Input.h>
#include <vsg/io/Output.h>

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

    void read(vsg::Input& input)
    {
        input.read("width", width);
        input.read("exp_term", exp_term);
        input.read("exp_scale", exp_scale);
        input.read("linear_term", linear_term);
        input.read("constant_term", constant_term);
    }

    void write(vsg::Output& output) const
    {
        output.write("width", width);
        output.write("exp_term", exp_term);
        output.write("exp_scale", exp_scale);
        output.write("linear_term", linear_term);
        output.write("constant_term", constant_term);
    }

};

class AtmosphereModelSettings : public vsg::Inherit<vsg::Object, AtmosphereModelSettings>
{
public:
    AtmosphereModelSettings(vsg::ref_ptr<vsg::EllipsoidModel> model);
    AtmosphereModelSettings();
    virtual ~AtmosphereModelSettings();

    void read(vsg::Input& input) override;
    void write(vsg::Output& output) const override;

    int numThreads = 8;

    int transmittanceWidth = 256;
    int transmittanceHeight = 64;

    int scaterringR = 32;
    int scaterringMU = 128;
    int scaterringMU_S = 32;
    int scaterringNU = 8;

    int irradianceWidth = 64;
    int irradianceHeight = 16;

    int precomputedWavelenghts = 15;

    float sunAngularRadius = 0.01935f;

    float miePhaseFunction_g = 0.8f;

    double maxSunZenithAngle = 120.0 / 180.0 * vsg::PI;
    double lengthUnitInMeters = 1000.0;

    vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel;

    /// Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
    /// (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
    /// summed and averaged in each bin (e.g. the value for 360nm is the average
    /// of the ASTM G-173 values for all wavelengths between 360 and 370nm).
    /// Values in W.m^-2.
    double kSolarIrradiance[48] = {
        1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
        1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
        1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
        1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
        1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
        1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
    };
    /// Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
    /// referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
    /// each bin (e.g. the value for 360nm is the average of the original values
    /// for all wavelengths between 360 and 370nm). Values in m^2.
    double kOzoneCrossSection[48] = {
        1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
        8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
        1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
        4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
        2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
        6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
        2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
    };
    /// From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
    double kDobsonUnit = 2.687e20;
    /// Maximum number density of ozone molecules, in m^-3 (computed so at to get
    /// 300 Dobson units of ozone - for this we divide 300 DU by the integral of
    /// the ozone density profile defined below, which is equal to 15km).
    double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
    /// Wavelength independent solar irradiance "spectrum" (not physically
    /// realistic, but was used in the original implementation).
    double kConstantSolarIrradiance = 1.5;

    double kAtmoshpereHeight = 60000.0;

    double kRayleigh = 1.24062e-6;
    double kRayleighScaleHeight = 8000.0;
    double kMieScaleHeight = 1200.0;
    double kMieAngstromAlpha = 0.0;
    double kMieAngstromBeta = 5.328e-3;
    double kMieSingleScatteringAlbedo = 0.9;
    double kGroundAlbedo = 0.1;

    bool use_constant_solar_spectrum_ = false;
    bool use_ozone = true;

    DensityProfileLayer rayleighDensityLayer = {0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0, lengthUnitInMeters};
    DensityProfileLayer mieDensityLayer = {0.0f, 1.0f, -1.0f / kMieScaleHeight, 0.0f, 0.0f, lengthUnitInMeters};

    DensityProfileLayer absorptionDensityLayer0 = {25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0, lengthUnitInMeters};
    DensityProfileLayer absorptionDensityLayer1 = {0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0, lengthUnitInMeters};
};

struct RuntimeSettings
{
    vsg::vec4 whitePoint;
    vsg::vec2 sunSize;

    void read(vsg::Input& input)
    {
        input.read("whitePoint", whitePoint);
        input.read("sunSize", sunSize);
    }

    void write(vsg::Output& output) const
    {
        output.write("whitePoint", whitePoint);
        output.write("sunSize", sunSize);
    }
};

class AtmosphereData : public vsg::Inherit<vsg::Object, AtmosphereData>
{
public:
    AtmosphereData();
    virtual ~AtmosphereData();

    void read(vsg::Input& input) override;
    void write(vsg::Output& output) const override;

    uint32_t cubeSize = 1024;
    int numViewerThreads = 32;

    double lengthUnitInMeters = 1000.0;

    vsg::ref_ptr<vsg::ImageInfo> transmittanceTexture;
    vsg::ref_ptr<vsg::ImageInfo> irradianceTexture;
    vsg::ref_ptr<vsg::ImageInfo> scatteringTexture;
    vsg::ref_ptr<vsg::ImageInfo> singleMieScatteringTexture;

    vsg::ref_ptr<vsg::Data> transmittanceData;
    vsg::ref_ptr<vsg::Data> irradianceData;
    vsg::ref_ptr<vsg::Data> scatteringData;
    vsg::ref_ptr<vsg::Data> singleMieScatteringData;

    vsg::ref_ptr<vsg::ImageInfo> skyBox;

    vsg::ref_ptr<vsg::ImageInfo> environmentMap;
    vsg::ref_ptr<vsg::ImageInfo> reflectionMap;

    vsg::ref_ptr<vsg::ShaderStage> reflectionMapShader;
    vsg::ref_ptr<vsg::ShaderStage> environmentMapShader;
    vsg::ref_ptr<vsg::ShaderSet> phongShaderSet;
    vsg::ref_ptr<vsg::ShaderSet> pbrShaderSet;
    vsg::ref_ptr<vsg::Node> sky;

    vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel;

    vsg::ref_ptr<vsg::Value<RuntimeSettings>> settings;

    vsg::dvec3 sunDirection = {0.0, std::sin(vsg::PI), std::cos(vsg::PI)};

    void setSunAngle(double radians);

    void setDate(tm time);

    vsg::ref_ptr<vsg::CommandGraph> createCubeMapGraph(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::vec4Value> camera);
    vsg::ref_ptr<vsg::View> createSkyView(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Camera> camera);

    vsg::ref_ptr<vsg::ImageInfo> createCubemap(uint32_t size);
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
    /// The distance between the bottom and the top of the atmosphere in m.
    /// </summary>
    double atmoshpereHeight;

    vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel;

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

    vsg::ref_ptr<vsg::DescriptorSetLayout> viewDescriptorSetLayout;

    vsg::ref_ptr<vsg::ShaderCompileSettings> compileSettings;

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


    vsg::ShaderStage::SpecializationConstants _computeConstants;
    vsg::ShaderStage::SpecializationConstants _renderConstants;

public:
    AtmosphereModel(vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Device> device, vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice, vsg::ref_ptr<vsg::Options> options);
    AtmosphereModel(vsg::ref_ptr<vsg::Device> device, vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice, vsg::ref_ptr<vsg::Options> options);
    virtual ~AtmosphereModel();

    void initialize(int scatteringOrders);

    vsg::vec3 convertSpectrumToLinearSrgb(double c);

    vsg::ref_ptr<AtmosphereData> getData();

private:
    void generateTextures();

    vsg::ref_ptr<vsg::ShaderSet> createPhongShaderSet();
    vsg::ref_ptr<vsg::Node> createSky();

    vsg::ref_ptr<vsg::ImageInfo> generate2D(uint32_t width, uint32_t height, bool init = false);
    vsg::ref_ptr<vsg::ImageInfo> generate3D(uint32_t width, uint32_t height, uint32_t depth, bool init = false);

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

    vsg::ref_ptr<vsg::Data> copyAndMapData(vsg::ref_ptr<vsg::ImageInfo> info);

    friend class AtmosphereLighting;
};

extern vsg::ref_ptr<AtmosphereModel> createAtmosphereModel(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Options> options);
extern vsg::ref_ptr<AtmosphereModel> createAtmosphereModel(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::EllipsoidModel> eps, vsg::ref_ptr<vsg::Options> options);
extern vsg::ref_ptr<AtmosphereModel> createAtmosphereModel(vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Options> options);
}

EVSG_type_name(atmosphere::AtmosphereModel)
EVSG_type_name(atmosphere::AtmosphereModelSettings)
EVSG_type_name(atmosphere::AtmosphereData)
EVSG_type_name(atmosphere::RuntimeSettings)

namespace vsg {
    template<>
    constexpr bool has_read_write<atmosphere::DensityProfileLayer>() { return true; }

    template<>
    constexpr bool has_read_write<atmosphere::RuntimeSettings>() { return true; }
}
#endif // ATMOSPHEREMODEL_H
