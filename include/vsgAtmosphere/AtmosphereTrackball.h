#ifndef ATMOSPHERETRACKBALL_H
#define ATMOSPHERETRACKBALL_H

#include <vsg/app/Trackball.h>

namespace atmosphere
{
    class AtmosphereTrackball : public vsg::Inherit<vsg::Trackball, AtmosphereTrackball>
    {
    public:
        explicit AtmosphereTrackball(vsg::ref_ptr<vsg::Camera> camera, vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel);

        void clampToGlobe() override;

        double minimumAltitude = 1.0;
        double maximumAltitude = 6000.0;
    };
}

#endif // ATMOSPHERETRACKBALL_H
