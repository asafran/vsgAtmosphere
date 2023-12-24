#include "AtmosphereTrackball.h"




atmosphere::AtmosphereTrackball::AtmosphereTrackball(vsg::ref_ptr<vsg::Camera> camera, vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel)
    : Inherit(camera, ellipsoidModel)
{

}

void atmosphere::AtmosphereTrackball::clampToGlobe()
{
    // get the location of the current lookAt center
    auto location_center = _ellipsoidModel->convertECEFToLatLongAltitude(_lookAt->center);
    auto location_eye = _ellipsoidModel->convertECEFToLatLongAltitude(_lookAt->eye);

    double ratio = location_eye.z / (location_eye.z - location_center.z);
    auto location = _ellipsoidModel->convertECEFToLatLongAltitude(_lookAt->center * ratio + _lookAt->eye * (1.0 - ratio));

    // clamp to the globe
    location.z = minimumAltitude;

    // compute clamped position back in ECEF
    auto ecef = _ellipsoidModel->convertLatLongAltitudeToECEF(location);

    // apply the new clamped position to the LookAt.
    _lookAt->center = ecef;

    if (location_eye.z < minimumAltitude)
    {
        location_eye.z = minimumAltitude;
        _lookAt->eye = _ellipsoidModel->convertLatLongAltitudeToECEF(location_eye);
        _thrown = false;
    } else if(location_eye.z > maximumAltitude)
    {
        location_eye.z = maximumAltitude;
        _lookAt->eye = _ellipsoidModel->convertLatLongAltitudeToECEF(location_eye);
        _thrown = false;
    }
}
