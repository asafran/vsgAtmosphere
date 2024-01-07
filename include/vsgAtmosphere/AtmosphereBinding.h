#ifndef ATMOSPHEREBINDING_H
#define ATMOSPHEREBINDING_H

#include "AtmoshpereConstants.h"
#include "AtmosphereImage.h"

#include <vsg/utils/ShaderSet.h>
#include <vsg/state/DescriptorSet.h>

namespace atmosphere {

class AtmosphereModelSettings;

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

    struct Positional
    {
        vsg::vec3 sunDirection;
        float exposure = 0.0f;
        vsg::vec3 cameraPos;
        float raduis = 0.0f;
        float mu_s = 0.0f;
    };

    class PositionalBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, PositionalBinding>
    {
    public:
        explicit PositionalBinding(uint32_t in_set = POSITIONAL_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;
        vsg::ref_ptr<vsg::DescriptorSet> descriptorSet;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;

        vsg::ref_ptr<vsg::Value<Positional>> positional;
    };

    class AtmosphereBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, AtmosphereBinding>
    {
    public:
        explicit AtmosphereBinding(uint32_t in_set = ATMOSHPERE_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;

        vsg::ref_ptr<Image> transmittanceTexture;
        vsg::ref_ptr<Image> irradianceTexture;
        vsg::ref_ptr<Image> scatteringTexture;
        vsg::ref_ptr<Image> singleMieScatteringTexture;

        bool radiance = false;

        vsg::ref_ptr<vsg::Value<RuntimeSettings>> settings;
    };

    struct CloudRuntimeSettings
    {
        float 	  cloudMinHeight = 4000.0f;
        float 	  cloudMaxHeight = 1500.0f;
        float 	  shapeNoiseScale = 0.00241f;
        float 	  detailNoiseScale = 0.01326f;

        float 	  detailNoiseModifier = 0.5f;
        float 	  turbulenceNoiseScale = 0.01793f;
        float 	  turbulenceAmount = 1.0f;
        float 	  cloudCoverage = 0.7f;

        vsg::vec3 	  windDirection = {1.0f, 0.0f, 0.0f};
        float	  windSpeed = 0.0f;

        float	  windShearOffset = 500.0f;
        uint16_t	  maxNumSteps = 8;
        float 	  lightStepLength = 128.0f;
        float 	  lightConeRadius = 0.4f;

        vsg::vec3      cloudBaseColor = {1.0f, 1.0f, 1.0f};
        float 	  precipitation = 1.0f;

        float 	  ambientLightFactor = 0.12f;
        float 	  sunLightFactor = 1.0f;
        float 	  henyeyGreensteinGForward = 0.179f;
        float 	  henyeyGreensteinGBackward = 0.6f;

        void read(vsg::Input& input)
        {
            /*input.read("density_to_sigma_s", density_to_sigma_s);
            input.read("phase_g", phase_g);
            input.read("density_to_sigma_t", density_to_sigma_t);
            input.read("primary_ray_marching_steps", primary_ray_marching_steps);
            input.read("box_min", box_min);
            input.read("secondary_ray_marching_steps", secondary_ray_marching_steps);
            input.read("box_max", box_max);
            input.read("enable_multi_scattering", enable_multi_scattering);
            input.read("g_c", g_c);
            input.read("g_d", g_d);
            input.read("wc0", wc0);
            input.read("wc1", wc1);
            input.read("wh", wh);
            input.read("shape_tile", shape_tile);
            input.read("detail_tile", detail_tile);
            input.read("blend_alpha", blend_alpha);
            input.read("cirrus", cirrus);
            input.read("cumulus", cumulus);*/
        }

        void write(vsg::Output& output) const
        {
            /*output.write("density_to_sigma_s", density_to_sigma_s);
            output.write("phase_g", phase_g);
            output.write("density_to_sigma_t", density_to_sigma_t);
            output.write("primary_ray_marching_steps", primary_ray_marching_steps);
            output.write("box_min", box_min);
            output.write("secondary_ray_marching_steps", secondary_ray_marching_steps);
            output.write("box_max", box_max);
            output.write("enable_multi_scattering", enable_multi_scattering);
            output.write("g_c", g_c);
            output.write("g_d", g_d);
            output.write("wc0", wc0);
            output.write("wc1", wc1);
            output.write("wh", wh);
            output.write("shape_tile", shape_tile);
            output.write("detail_tile", detail_tile);
            output.write("blend_alpha", blend_alpha);
            output.write("cirrus", cirrus);
            output.write("cumulus", cumulus);*/
        }
    };

    struct Parameters
    {
        vsg::vec4 solar_irradiance;
        vsg::vec4 rayleigh_scattering;
        vsg::vec4 mie_scattering;
        vsg::vec4 mie_extinction;
        vsg::vec4 ground_albedo;
        vsg::vec4 absorption_extinction;
    };

    class ComputeParametersBinding : public vsg::Inherit<vsg::Object, ComputeParametersBinding>
    {
    public:
        explicit ComputeParametersBinding(vsg::ref_ptr<AtmosphereModelSettings> settings);

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;
        vsg::ref_ptr<vsg::DescriptorSet> descriptorSet;

        void computeLuminanceFromRadiance(const vsg::vec3 &lambdas, double dlambda);
        void resetLuminanceFromRadiance();
        void setParameters(const Parameters& value);

    private:

        vsg::ref_ptr<vsg::mat4Value> _luminacneFromRadiance;
        vsg::ref_ptr<vsg::DescriptorBuffer> _luminacneFromRadianceBuffer;
        vsg::ref_ptr<vsg::Value<Parameters>> _parameters;
        vsg::ref_ptr<vsg::DescriptorBuffer> _parametersBuffer;
        vsg::ref_ptr<vsg::Array<DensityProfileLayer>> _profileLayers;
        vsg::ref_ptr<vsg::DescriptorBuffer> _profileLayersBuffer;

    };

    class CloudsBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, CloudsBinding>
    {
    public:
        explicit CloudsBinding(uint32_t in_set = CLOUDS_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;

        vsg::ref_ptr<Image> shapeNoiseTexture;
        vsg::ref_ptr<Image> detailNoiseTexture;
        vsg::ref_ptr<Image> blueNoiseTexture;
        vsg::ref_ptr<Image> curlNoiseTexture;

        vsg::ref_ptr<vsg::Value<CloudRuntimeSettings>> settings;
    };

    class BRDFBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, BRDFBinding>
    {
    public:
        explicit BRDFBinding(uint32_t in_set = PBR_DESCRIPTOR_SET);
        BRDFBinding(unsigned int size, unsigned int layers, uint32_t in_set = PBR_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;

        vsg::ref_ptr<Image> BRDFlutTexture;
        vsg::ref_ptr<vsg::ImageInfo> prefilteredTexture;

    private:
        void createCubemap(uint32_t size, uint32_t layers);
    };
}

#endif // ATMOSPHEREBINDING_H
