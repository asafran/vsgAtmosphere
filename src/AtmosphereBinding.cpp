#include "AtmosphereBinding.h"
#include "Atmosphere.h"
#include "AtmosphereTools.h"
#include "vsg/state/BindDescriptorSet.h"
#include "vsg/state/Descriptor.h"
#include "vsg/state/DescriptorBuffer.h"
#include "vsg/state/DescriptorImage.h"
#include "vsg/state/DescriptorSet.h"

namespace atmosphere {

    PositionalBinding::PositionalBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        vsg::DescriptorSetLayoutBindings descriptorBindings{{0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}};
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        positional = vsg::Value<Positional>::create();
        positional->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;

        descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, vsg::Descriptors{{vsg::DescriptorBuffer::create(positional)}});
    }

    int PositionalBinding::compare(const Object &rhs) const
    {
        return CustomDescriptorSetBinding::compare(rhs);
    }

    bool PositionalBinding::compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout &dsl) const
    {
        return descriptorSetLayout->compare(dsl) == 0;
    }

    vsg::ref_ptr<vsg::DescriptorSetLayout> PositionalBinding::createDescriptorSetLayout()
    {
        return descriptorSetLayout;
    }

    vsg::ref_ptr<vsg::StateCommand> PositionalBinding::createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout)
    {
        return vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_GRAPHICS, layout, set, descriptorSet);
    }

    AtmosphereBinding::AtmosphereBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        vsg::DescriptorSetLayoutBindings descriptorBindings{
           {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {5, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}};
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        settings = vsg::Value<RuntimeSettings>::create();
        settings->properties.dataVariance = vsg::DYNAMIC_DATA;
    }

    int AtmosphereBinding::compare(const Object &rhs) const
    {
        return CustomDescriptorSetBinding::compare(rhs);
    }

    void AtmosphereBinding::read(vsg::Input &input)
    {
        input.readObject("transmittanceData", transmittanceTexture);
        input.readObject("irradianceData", irradianceTexture);
        input.readObject("scatteringData", scatteringTexture);
        input.readObject("singleMieScatteringData", singleMieScatteringTexture);
    }

    void AtmosphereBinding::write(vsg::Output &output) const
    {
        output.writeObject("transmittanceData", transmittanceTexture);
        output.writeObject("irradianceData", irradianceTexture);
        output.writeObject("scatteringData", scatteringTexture);
        output.writeObject("singleMieScatteringData", singleMieScatteringTexture);
    }

    bool AtmosphereBinding::compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout &dsl) const
    {
        return descriptorSetLayout->compare(dsl) == 0;
    }

    vsg::ref_ptr<vsg::DescriptorSetLayout> AtmosphereBinding::createDescriptorSetLayout()
    {
        return descriptorSetLayout;
    }

    vsg::ref_ptr<vsg::StateCommand> AtmosphereBinding::createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout)
    {
        vsg::Descriptors descriptors{
            {vsg::DescriptorImage::create(transmittanceTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(irradianceTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(scatteringTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(singleMieScatteringTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorBuffer::create(settings, 4, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER)}
        };
        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        auto bindDescriptorSet = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_GRAPHICS, layout, set, descriptorSet);
        return bindDescriptorSet;
    }

    CloudsBinding::CloudsBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}
        };
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        settings = vsg::Value<CloudRuntimeSettings>::create();
        settings->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;
    }

    int CloudsBinding::compare(const Object &rhs) const
    {
        return CustomDescriptorSetBinding::compare(rhs);
    }

    void CloudsBinding::read(vsg::Input &input)
    {
        input.readObject("shapeNoiseTexture", shapeNoiseTexture);
        input.readObject("detailNoiseTexture", detailNoiseTexture);
        input.readObject("blueNoiseTexture", blueNoiseTexture);
        input.readObject("curlNoiseTexture", curlNoiseTexture);
    }

    void CloudsBinding::write(vsg::Output &output) const
    {
        output.writeObject("shapeNoiseTexture", shapeNoiseTexture);
        output.writeObject("detailNoiseTexture", detailNoiseTexture);
        output.writeObject("blueNoiseTexture", blueNoiseTexture);
        output.writeObject("curlNoiseTexture", curlNoiseTexture);
    }

    bool CloudsBinding::compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout &dsl) const
    {
        return descriptorSetLayout->compare(dsl) == 0;
    }

    vsg::ref_ptr<vsg::DescriptorSetLayout> CloudsBinding::createDescriptorSetLayout()
    {
        return descriptorSetLayout;
    }

    vsg::ref_ptr<vsg::StateCommand> CloudsBinding::createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout)
    {
        vsg::Descriptors descriptors{
            {vsg::DescriptorImage::create(shapeNoiseTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(detailNoiseTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(blueNoiseTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(curlNoiseTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorBuffer::create(settings, 4, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER)}
        };
        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        auto bindDescriptorSet = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_GRAPHICS, layout, set, descriptorSet);
        return bindDescriptorSet;
    }

    BRDFBinding::BRDFBinding(unsigned int size, unsigned int layers, uint32_t in_set)
        : Inherit(in_set)
    {
        createCubemap(size, layers);

        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}};
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);
    }

    void BRDFBinding::createCubemap(uint32_t size, uint32_t layers)
    {
        vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
        image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT);
        image->flags |= VK_IMAGE_CREATE_CUBE_COMPATIBLE_BIT;
        image->format = VK_FORMAT_R16G16B16A16_SFLOAT;
        image->mipLevels = layers;
        image->extent = VkExtent3D{size, size, 1};
        image->imageType = VK_IMAGE_TYPE_2D;
        image->arrayLayers = 6;

        auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
        auto sampler = vsg::Sampler::create();
        sampler->maxLod = layers;
        prefilteredTexture = vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
    }

    BRDFBinding::BRDFBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}};
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);
    }


    int BRDFBinding::compare(const Object &rhs) const
    {
        return CustomDescriptorSetBinding::compare(rhs);
    }


    void BRDFBinding::read(vsg::Input &input)
    {
        unsigned int size = 128;
        unsigned int layers = 6;

        input.readObject("BRDFlutTexture", BRDFlutTexture);

        createCubemap(size, layers);
    }


    void BRDFBinding::write(vsg::Output &output) const
    {
        output.writeObject("BRDFlutTexture", BRDFlutTexture);
    }


    bool BRDFBinding::compatibleDescriptorSetLayout(const vsg::DescriptorSetLayout &dsl) const
    {
        return descriptorSetLayout->compare(dsl) == 0;
    }


    vsg::ref_ptr<vsg::DescriptorSetLayout> BRDFBinding::createDescriptorSetLayout()
    {
        return descriptorSetLayout;
    }


    vsg::ref_ptr<vsg::StateCommand> BRDFBinding::createStateCommand(vsg::ref_ptr<vsg::PipelineLayout> layout)
    {
        vsg::Descriptors descriptors{
            {vsg::DescriptorImage::create(BRDFlutTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
            {vsg::DescriptorImage::create(prefilteredTexture, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)}
        };
        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        auto bindDescriptorSet = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_GRAPHICS, layout, set, descriptorSet);
        return bindDescriptorSet;
    }
    
    ComputeParametersBinding::ComputeParametersBinding(vsg::ref_ptr<AtmosphereModelSettings> settings)
    {
        _luminacneFromRadiance = vsg::mat4Value::create(vsg::mat4());
        _parameters = vsg::Value<Parameters>::create(Parameters());
        _profileLayers = vsg::Array<DensityProfileLayer>::create(4);
        _profileLayers->set(0, settings->rayleighDensityLayer);
        _profileLayers->set(1, settings->mieDensityLayer);
        _profileLayers->set(2, settings->absorptionDensityLayer0);
        _profileLayers->set(3, settings->absorptionDensityLayer1);

        _luminacneFromRadianceBuffer = vsg::DescriptorBuffer::create(_luminacneFromRadiance, 0, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER);
        _parametersBuffer = vsg::DescriptorBuffer::create(_parameters, 1, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER);
        _profileLayersBuffer = vsg::DescriptorBuffer::create(_profileLayers, 2, 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER);

        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
            {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        };
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        vsg::Descriptors descriptors{_luminacneFromRadianceBuffer, _parametersBuffer, _profileLayersBuffer};
        descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
    }

    void ComputeParametersBinding::computeLuminanceFromRadiance(const vsg::vec3 &lambdas, double dlambda)
    {

        auto coeff = [dlambda](double lambda, int component)
        {
            // Note that we don't include MAX_LUMINOUS_EFFICACY here, to avoid
            // artefacts due to too large values when using half precision on GPU.
            // We add this term back in kAtmosphereShader, via
            // SKY_SPECTRAL_RADIANCE_TO_LUMINANCE (see also the comments in the
            // Model constructor).
            double x = cieColorMatchingFunctionTableValue(lambda, 1);
            double y = cieColorMatchingFunctionTableValue(lambda, 2);
            double z = cieColorMatchingFunctionTableValue(lambda, 3);
            return static_cast<float>((
                                          XYZ_TO_SRGB[component * 3] * x +
                                          XYZ_TO_SRGB[component * 3 + 1] * y +
                                          XYZ_TO_SRGB[component * 3 + 2] * z) * dlambda);
        };
        vsg::mat4 luminance_from_radiance{
                                          coeff(lambdas[0], 0), coeff(lambdas[0], 1), coeff(lambdas[0], 2), 0.0f,
                                          coeff(lambdas[1], 0), coeff(lambdas[1], 1), coeff(lambdas[1], 2), 0.0f,
                                          coeff(lambdas[2], 0), coeff(lambdas[2], 1), coeff(lambdas[2], 2), 0.0f,
                                          0.0, 0.0f, 0.0f, 1.0f};

        _luminacneFromRadiance->set(luminance_from_radiance);
        _luminacneFromRadianceBuffer->copyDataListToBuffers();
    }

    void ComputeParametersBinding::resetLuminanceFromRadiance()
    {
        _luminacneFromRadiance->set(vsg::mat4());
        _luminacneFromRadianceBuffer->copyDataListToBuffers();
    }

    void ComputeParametersBinding::setParameters(const Parameters &value)
    {
        _parameters->set(value);
        _parametersBuffer->copyDataListToBuffers();
    }
}
