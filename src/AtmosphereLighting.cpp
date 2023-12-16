#include "AtmosphereLighting.h"
#include <vsg/app/View.h>
#include <vsg/io/read.h>
#include <vsg/state/material.h>
#include <vsg/state/DescriptorImage.h>

namespace atmosphere {

    AtmosphereLighting::AtmosphereLighting(vsg::ref_ptr<vsg::View> view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime)
        : Inherit(view.get())
    , atmosphereRuntime(in_atmosphereRuntime)
    {
        atmosphereRuntime->atmosphereBinding->createDescriptorSet();
    }

    void AtmosphereLighting::bindDescriptorSets(vsg::CommandBuffer &commandBuffer, VkPipelineBindPoint pipelineBindPoint, VkPipelineLayout layout, uint32_t firstSet)
    {
        vsg::ViewDependentState::bindDescriptorSets(commandBuffer, pipelineBindPoint, layout, firstSet);
        
        auto vk = atmosphereRuntime->atmosphereBinding->descriptorSet->vk(commandBuffer.deviceID);
        vkCmdBindDescriptorSets(commandBuffer, pipelineBindPoint, layout, atmosphereRuntime->atmosphereBinding->set, 1, &vk, 0, nullptr);
    }

    void AtmosphereLighting::compile(vsg::Context &context)
    {
        vsg::ViewDependentState::compile(context);

        atmosphereRuntime->atmosphereBinding->descriptorSet->compile(context);
    }

    void AtmosphereLighting::traverse(vsg::RecordTraversal& rt) const
    {
        atmosphereRuntime->atmosphereBinding->positional->dirty();
        Positional p;
        auto mv =  view->camera->viewMatrix->transform();
        auto eye_direction = -normalize(atmosphereRuntime->sunDirection * vsg::inverse_3x3(mv));
        vsg::vec4 directionExposure{static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z), static_cast<float>(exposure * 1e-6)};
        p.sunDirectionExp = directionExposure;

        auto eye_position = mv * vsg::dvec3();
        eye_position /= atmosphereRuntime->lengthUnitInMeters;
        p.cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z), 0.0f);
        atmosphereRuntime->atmosphereBinding->positional->set(p);

        vsg::ViewDependentState::traverse(rt);
    }

    AtmosphereLighting::~AtmosphereLighting()
    {
    }

    SkyLighting::SkyLighting(vsg::ref_ptr<vsg::View> view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime)
        : Inherit(view, in_atmosphereRuntime)
    {
        if(atmosphereRuntime->cloudsBinding)
            atmosphereRuntime->cloudsBinding->createDescriptorSet();
    }


    void SkyLighting::bindDescriptorSets(vsg::CommandBuffer &commandBuffer, VkPipelineBindPoint pipelineBindPoint, VkPipelineLayout layout, uint32_t firstSet)
    {
        AtmosphereLighting::bindDescriptorSets(commandBuffer, pipelineBindPoint, layout, firstSet);

        if(atmosphereRuntime->cloudsBinding)
        {
            auto vk = atmosphereRuntime->cloudsBinding->descriptorSet->vk(commandBuffer.deviceID);
            vkCmdBindDescriptorSets(commandBuffer, pipelineBindPoint, layout, atmosphereRuntime->cloudsBinding->set, 1, &vk, 0, nullptr);
        }
    }


    void SkyLighting::compile(vsg::Context &context)
    {
        AtmosphereLighting::compile(context);
        if(atmosphereRuntime->cloudsBinding)
            atmosphereRuntime->cloudsBinding->descriptorSet->compile(context);
    }


    void SkyLighting::traverse(vsg::RecordTraversal &rt) const
    {
        Positional p;
        auto mv =  view->camera->viewMatrix->transform();
        auto eye_direction = -normalize(-atmosphereRuntime->sunDirection);
        vsg::vec4 directionExposure{static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z), static_cast<float>(exposure * 1e-6)};
        p.sunDirectionExp = directionExposure;


        auto eye_position = mv[3].xyz;
        eye_position /= atmosphereRuntime->lengthUnitInMeters;
        p.cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z), 0.0f);
        //atmosphereRuntime->atmosphereBinding->positional->set(p);
        //atmosphereRuntime->atmosphereBinding->positional->dirty();

        vsg::ViewDependentState::traverse(rt);
    }


    SkyLighting::~SkyLighting()
    {

    }

}
