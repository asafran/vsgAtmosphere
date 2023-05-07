#include "AtmosphereLighting.h"
#include <vsg/io/read.h>
#include <vsg/state/material.h>
#include <vsg/state/DescriptorImage.h>

namespace atmosphere {

    AtmosphereLighting::AtmosphereLighting(AtmosphereModel *model, vsg::ref_ptr<vsg::ViewMatrix> view, uint32_t maxNumberLights, uint32_t maxViewports)
        : vsg::Inherit<vsg::ViewDependentState, AtmosphereLighting>(maxNumberLights, maxViewports)
        , atmosphereModel(model)
        , viewMatrix(view)
        , positional(vsg::Value<Positional>::create())
    {
        positional->properties.dataVariance = vsg::DYNAMIC_DATA_TRANSFER_AFTER_RECORD;

        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(model->_transmittanceTexture, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(model->_irradianceTexture, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(model->_scatteringTexture, 4, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorImage::create(model->_singleMieScatteringTexture, 5, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
        descriptorSet->descriptors.push_back(vsg::DescriptorBuffer::create(model->runtimeSettings, 6, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER));
        descriptorSet->descriptors.push_back(vsg::DescriptorBuffer::create(positional, 7, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER));

        descriptorSetLayout->bindings.push_back({2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({5, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({6, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
        descriptorSetLayout->bindings.push_back({7, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, nullptr});
    }

    void AtmosphereLighting::pack()
    {
        auto mv = viewMatrix.valid() ? viewMatrix->transform() : vsg::dmat4();
        auto eye_direction = -normalize(atmosphereModel->sunDirection * vsg::inverse_3x3(mv));
        positional->value().sunDirection.set(static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z), 0.0f);

        auto eye_global = normalize(atmosphereModel->sunDirection);
        positional->value().globalSunDirection.set(static_cast<float>(eye_global.x), static_cast<float>(eye_global.y), static_cast<float>(eye_global.z), 0.0f);

        auto eye_position =  mv * vsg::dvec3();
        eye_position /= 1000.0;
        positional->value().cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z), 0.0f);
        positional->dirty();

        vsg::ViewDependentState::pack();
    }

    AtmosphereLighting::~AtmosphereLighting()
    {
    }
}
