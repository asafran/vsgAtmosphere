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
        AtmosphereLighting(vsg::ref_ptr<vsg::View> in_view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime);

        vsg::ref_ptr<AtmosphereRuntime> atmosphereRuntime;

        double exposure = 10.0;

        void traverse(vsg::Visitor& visitor) override { ViewDependentState::traverse(visitor); atmosphereRuntime->accept(visitor); }
        void traverse(vsg::ConstVisitor& visitor) const override { ViewDependentState::traverse(visitor); atmosphereRuntime->accept(visitor); }
        void traverse(vsg::RecordTraversal& rt) const override;

    protected:
        ~AtmosphereLighting();
    };

    class SkyLighting : public vsg::Inherit<AtmosphereLighting, SkyLighting>
    {
    public:
        SkyLighting(vsg::ref_ptr<vsg::View> in_view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime);

        void traverse(vsg::RecordTraversal& rt) const override;

    protected:
        ~SkyLighting();
    };
}

#endif // ATMOSPHERELIGHTING_H
