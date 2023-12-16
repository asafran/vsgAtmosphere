#include "AtmosphereBinding.h"
#include "vsg/state/Descriptor.h"
#include "vsg/state/DescriptorBuffer.h"
#include "vsg/state/DescriptorImage.h"
#include "vsg/state/DescriptorSet.h"

namespace atmosphere {

    AtmosphereBinding::AtmosphereBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        vsg::DescriptorSetLayoutBindings descriptorBindings{
           {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
           {5, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}};
        descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        settings = vsg::Value<RuntimeSettings>::create();
        settings->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;
        positional = vsg::Value<Positional>::create(Positional{{1.0, 0.0, 1.0, 0.0},{0.0, 1.0, 0.0, 1.0}});
        positional->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;
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
        return {};
    }

    vsg::ref_ptr<vsg::DescriptorSet> AtmosphereBinding::createDescriptorSet()
    {
        if(!descriptorSet)
        {
            vsg::Descriptors descriptors{
                {vsg::DescriptorImage::create(transmittanceTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(irradianceTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(scatteringTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(singleMieScatteringTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorBuffer::create(settings, 4, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER)},
                {vsg::DescriptorBuffer::create(positional, 5, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER)},
            };
            descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        }
        return descriptorSet;
    }

    CloudsBinding::CloudsBinding(uint32_t in_set)
        : Inherit(in_set)
    {
        // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr},
            {4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr}
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
        return {};
    }

    vsg::ref_ptr<vsg::DescriptorSet> CloudsBinding::createDescriptorSet()
    {
        if(!descriptorSet)
        {
            vsg::Descriptors descriptors{
                {vsg::DescriptorImage::create(shapeNoiseTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(detailNoiseTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(blueNoiseTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorImage::create(curlNoiseTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER)},
                {vsg::DescriptorBuffer::create(settings, 4, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER)}
            };
            descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        }
        return descriptorSet;
    }
}
