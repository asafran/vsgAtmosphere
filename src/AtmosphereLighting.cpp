#include "AtmosphereLighting.h"
#include <vsg/io/read.h>
#include <vsg/state/material.h>
#include <vsg/state/DescriptorImage.h>

namespace atmosphere {

    AtmosphereLighting::AtmosphereLighting(vsg::ref_ptr<vsg::ViewMatrix> view, uint32_t maxNumberLights, uint32_t maxViewports)
        : vsg::Inherit<vsg::ViewDependentState, AtmosphereLighting>(maxNumberLights, maxViewports)
        , _viewMatrix(view)
        , _positional(vsg::Value<Positional>::create())
    {
        _positional->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;

        descriptorSet->descriptors.push_back(vsg::DescriptorBuffer::create(_positional, 7, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER));

        descriptorSetLayout->bindings.push_back({2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({5, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({6, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({7, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
    }

    void AtmosphereLighting::assignData(vsg::ref_ptr<AtmosphereData> data)
    {
        if(_atmosphereData)
            return;
        _atmosphereData = data;

        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(data->transmittanceTexture, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(data->irradianceTexture, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(data->scatteringTexture, 4, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(data->singleMieScatteringTexture, 5, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorBuffer::create(data->settings, 6, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER));
    }

    void AtmosphereLighting::pack()
    {
        auto mv = _viewMatrix->transform();
        auto eye_direction = -normalize(transform ? _atmosphereData->sunDirection * vsg::inverse_3x3(mv) : -_atmosphereData->sunDirection);
        vsg::vec4 directionExposure{static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z), static_cast<float>(exposure * 1e-6)};
        _positional->value().sunDirectionExp = directionExposure;


        auto eye_position = transform ? mv * vsg::dvec3() : mv[3].xyz;
        eye_position /= _atmosphereData->lengthUnitInMeters;
        _positional->value().cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z), 0.0f);
        _positional->dirty();

        vsg::ViewDependentState::pack();
    }

    AtmosphereLighting::~AtmosphereLighting()
    {
    }
}
