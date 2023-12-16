#ifndef ATMOSPHERELIGHTING_H
#define ATMOSPHERELIGHTING_H

#include "Atmosphere.h"
#include "AtmosphereRuntime.h"
#include <vsg/state/ViewDependentState.h>
#include <vsg/utils/ShaderSet.h>

namespace atmosphere {

    class AtmosphereLighting : public vsg::Inherit<vsg::ViewDependentState, AtmosphereLighting>
    {
    public:
        AtmosphereLighting(vsg::ref_ptr<vsg::View> view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime);

        vsg::ref_ptr<AtmosphereRuntime> atmosphereRuntime;

        double exposure = 10.0;

        void bindDescriptorSets(vsg::CommandBuffer& commandBuffer, VkPipelineBindPoint pipelineBindPoint, VkPipelineLayout layout, uint32_t firstSet) override;
        void compile(vsg::Context& context) override;

        void traverse(vsg::RecordTraversal& rt) const override;

    protected:
        ~AtmosphereLighting();
    };

    class SkyLighting : public vsg::Inherit<AtmosphereLighting, SkyLighting>
    {
    public:
        SkyLighting(vsg::ref_ptr<vsg::View> view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime);

        void bindDescriptorSets(vsg::CommandBuffer& commandBuffer, VkPipelineBindPoint pipelineBindPoint, VkPipelineLayout layout, uint32_t firstSet) override;
        void compile(vsg::Context& context) override;

        void traverse(vsg::RecordTraversal& rt) const override;

    protected:
        ~SkyLighting();
    };
}

#endif // ATMOSPHERELIGHTING_H
