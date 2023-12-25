#include "AtmosphereLighting.h"
#include <vsg/app/View.h>
#include <vsg/io/read.h>
#include <vsg/state/material.h>
#include <vsg/state/DescriptorImage.h>

namespace atmosphere {

    AtmosphereLighting::AtmosphereLighting(vsg::ref_ptr<vsg::View> in_view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime)
        : Inherit(in_view.get())
        , atmosphereRuntime(in_atmosphereRuntime)
    {
    }

    void AtmosphereLighting::traverse(vsg::RecordTraversal& rt) const
    {

        Positional p;
        auto mv =  view->camera->viewMatrix->transform();
        auto eye_direction = -normalize(atmosphereRuntime->sunDirection * vsg::inverse_3x3(mv));
        p.sunDirection.set(static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z));
        p.exposure = static_cast<float>(exposure * atmosphereRuntime->exposureModifier);

        auto eye_position = mv * vsg::dvec3();
        eye_position /= atmosphereRuntime->lengthUnitInMeters;
        p.cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z));
        p.raduis = static_cast<float>(vsg::length(eye_position));

        p.mu_s = vsg::dot(eye_position, eye_direction) / p.raduis;
        atmosphereRuntime->positionalBinding->positional->set(p);
        atmosphereRuntime->positionalBinding->positional->dirty();

        vsg::ViewDependentState::traverse(rt);
    }

    AtmosphereLighting::~AtmosphereLighting()
    {
    }

    SkyLighting::SkyLighting(vsg::ref_ptr<vsg::View> in_view, vsg::ref_ptr<AtmosphereRuntime> in_atmosphereRuntime)
        : Inherit(in_view, in_atmosphereRuntime)
    {
    }

    void SkyLighting::traverse(vsg::RecordTraversal &) const
    {
        Positional p;
        auto mv =  view->camera->viewMatrix->transform();
        auto eye_direction = -normalize(-atmosphereRuntime->sunDirection);
        p.sunDirection.set(static_cast<float>(eye_direction.x), static_cast<float>(eye_direction.y), static_cast<float>(eye_direction.z));
        p.exposure = static_cast<float>(exposure * atmosphereRuntime->exposureModifier);

        auto eye_position = mv[3].xyz;
        eye_position /= atmosphereRuntime->lengthUnitInMeters;
        p.cameraPos.set(static_cast<float>(eye_position.x), static_cast<float>(eye_position.y), static_cast<float>(eye_position.z));
        p.raduis = static_cast<float>(vsg::length(eye_position));

        p.mu_s = vsg::dot(eye_position, eye_direction) / p.raduis;
        atmosphereRuntime->inversePositionalBinding->positional->set(p);
        atmosphereRuntime->inversePositionalBinding->positional->dirty();
    }


    SkyLighting::~SkyLighting()
    {

    }

}
