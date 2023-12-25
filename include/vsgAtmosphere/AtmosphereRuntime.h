#ifndef ATMOSPHERERUNTIME_H
#define ATMOSPHERERUNTIME_H

#include "AtmosphereBinding.h"

#include <vsg/app/Window.h>
#include <vsg/maths/common.h>
#include <vsg/app/EllipsoidModel.h>

namespace atmosphere {

    class AtmosphereRuntime : public vsg::Inherit<vsg::Object, AtmosphereRuntime>
    {
    public:
        explicit AtmosphereRuntime(vsg::ref_ptr<AtmosphereBinding> atmosphere, vsg::ref_ptr<CloudsBinding> clouds = {});
        AtmosphereRuntime();

        int compare(const Object& rhs) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

        bool createPhongShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts, bool radiance);
        bool createPBRShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts, bool radiance);
        bool createSkyShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts, bool radiance);

        uint32_t cubeSize = 1024;
        int numViewerThreads = 32;

        double exposureModifier = 1e-6;

        double lengthUnitInMeters = 1000.0;

        vsg::ref_ptr<vsg::ShaderSet> skyShaderSet;
        vsg::ref_ptr<vsg::ShaderSet> phongShaderSet;
        vsg::ref_ptr<vsg::ShaderSet> pbrShaderSet;

        vsg::ref_ptr<vsg::EllipsoidModel> ellipsoidModel;

        vsg::ref_ptr<AtmosphereBinding> atmosphereBinding;
        vsg::ref_ptr<CloudsBinding> cloudsBinding;

        vsg::ref_ptr<PositionalBinding> positionalBinding;
        vsg::ref_ptr<PositionalBinding> inversePositionalBinding;

        vsg::dvec3 sunDirection = {0.0, std::sin(vsg::PI), std::cos(vsg::PI)};

        void setSunAngle(double radians);

        void setDate(tm time);

        vsg::ref_ptr<vsg::CommandGraph> createCubeMapGraph(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::vec4Value> camera);
        vsg::ref_ptr<vsg::View> createSkyView(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Camera> camera);
        vsg::ref_ptr<vsg::Node> createSky(bool viewerInSpace = false);
    };

}

#endif // ATMOSPHERERUNTIME_H
