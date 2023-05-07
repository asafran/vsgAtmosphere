#ifndef ATMOSPHERELIGHTING_H
#define ATMOSPHERELIGHTING_H

#include "Atmosphere.h"
#include <vsg/state/ViewDependentState.h>
#include <vsg/utils/ShaderSet.h>

namespace atmosphere {

    class AtmosphereLighting : public vsg::Inherit<vsg::ViewDependentState, AtmosphereLighting>
    {
    public:
        struct Positional
        {
            vsg::vec4 sunDirection;
            vsg::vec4 globalSunDirection;
            vsg::vec4 cameraPos;
        };

        AtmosphereLighting(AtmosphereModel *model, vsg::ref_ptr<vsg::ViewMatrix> view = {}, uint32_t maxNumberLights = 64, uint32_t maxViewports = 1);

        vsg::ref_ptr<AtmosphereModel> atmosphereModel;
        vsg::ref_ptr<vsg::ViewMatrix> viewMatrix;

        vsg::ref_ptr<vsg::Value<Positional>> positional;

        void pack() override;

    protected:
        ~AtmosphereLighting();
    };
}

#endif // ATMOSPHERELIGHTING_H
