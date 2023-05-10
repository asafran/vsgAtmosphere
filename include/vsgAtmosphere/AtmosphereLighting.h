#ifndef ATMOSPHERELIGHTING_H
#define ATMOSPHERELIGHTING_H

#include "Atmosphere.h"
#include <vsg/state/ViewDependentState.h>
#include <vsg/utils/ShaderSet.h>

namespace atmosphere {

    class AtmosphereLighting : public vsg::Inherit<vsg::ViewDependentState, AtmosphereLighting>
    {
    public:
        AtmosphereLighting(vsg::ref_ptr<vsg::ViewMatrix> view = {}, uint32_t maxNumberLights = 64, uint32_t maxViewports = 1);

        void assignData(vsg::ref_ptr<AtmosphereData> data);

        double exposure = 10.0;
        bool transform = true;

        void pack() override;

    protected:
        ~AtmosphereLighting();

        struct Positional
        {
            vsg::vec4 sunDirectionExp;
            vsg::vec4 cameraPos;
        };

        vsg::ref_ptr<AtmosphereData> _atmosphereData;
        vsg::ref_ptr<vsg::ViewMatrix> _viewMatrix;

        vsg::ref_ptr<vsg::Value<Positional>> _positional;
    };
}

#endif // ATMOSPHERELIGHTING_H
