#ifndef ATMOSPHEREBINDING_H
#define ATMOSPHEREBINDING_H

#include "AtmoshpereConstatnts.h"
#include "AtmosphereImage.h"

#include <vsg/utils/ShaderSet.h>

namespace atmosphere {

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
        vsg::vec4 sunDirectionExp;
        vsg::vec4 cameraPos;
    };

    class AtmosphereBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, AtmosphereBinding>
    {
    public:
        AtmosphereBinding(uint32_t in_set = ATMOSHPERE_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;
        vsg::ref_ptr<vsg::DescriptorSet> descriptorSet;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;
        vsg::ref_ptr<vsg::DescriptorSet> createDescriptorSet();

        vsg::ref_ptr<Image> transmittanceTexture;
        vsg::ref_ptr<Image> irradianceTexture;
        vsg::ref_ptr<Image> scatteringTexture;
        vsg::ref_ptr<Image> singleMieScatteringTexture;

        vsg::ref_ptr<vsg::Value<RuntimeSettings>> settings;
        vsg::ref_ptr<vsg::Value<Positional>> positional;
    };

    struct CloudRuntimeSettings
    {
        float 	  cloudMinHeight = 1500.0f;
        float 	  cloudMaxHeight = 4000.0f;
        float 	  shapeNoiseScale = 0.3f;
        float 	  detailNoiseScale = 5.5f;

        float 	  detailNoiseModifier = 0.5f;
        float 	  turbulenceNoiseScale = 7.440f;
        float 	  turbulenceAmount = 1.0f;
        float 	  cloudCoverage = 0.7f;

        vsg::vec3 	  windDirection = {1.0f, 0.0f, 0.0f};
        float	  windSpeed = 0.0f;

        float	  windShearOffset = 500.0f;
        uint16_t	  maxNumSteps = 128;
        float 	  lightStepLength = 64.0f;
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

    class CloudsBinding : public vsg::Inherit<vsg::CustomDescriptorSetBinding, CloudsBinding>
    {
    public:
        CloudsBinding(uint32_t in_set = CLOUDS_DESCRIPTOR_SET);

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout;
        vsg::ref_ptr<vsg::DescriptorSet> descriptorSet;

        bool compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout& dsl) const override;
        vsg::ref_ptr<vsg::DescriptorSetLayout> createDescriptorSetLayout() override;

        vsg::ref_ptr<vsg::StateCommand> createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout) override;
        vsg::ref_ptr<vsg::DescriptorSet> createDescriptorSet();

        vsg::ref_ptr<Image> shapeNoiseTexture;
        vsg::ref_ptr<Image> detailNoiseTexture;
        vsg::ref_ptr<Image> blueNoiseTexture;
        vsg::ref_ptr<Image> curlNoiseTexture;

        vsg::ref_ptr<vsg::Value<CloudRuntimeSettings>> settings;
    };
}

#endif // ATMOSPHEREBINDING_H
