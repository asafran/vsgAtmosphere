#ifndef CLOUDS_H
#define CLOUDS_H

#include "Atmosphere.h"
#include "AtmosphereBinding.h"
#include "AtmosphereImage.h"
#include <vsg/app/Window.h>
#include <vsg/vk/Device.h>
#include <vsg/state/ImageInfo.h>
#include <vsg/io/Options.h>

namespace atmosphere {


    class CloudsGenerator : public vsg::Inherit<vsg::Object, CloudsGenerator>
    {
    public:
        CloudsGenerator(vsg::ref_ptr<vsg::Device> device,
                        vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice,
                        vsg::ref_ptr<AtmosphereModelSettings> settings,
                        vsg::ref_ptr<vsg::Options> options);

        void initialize();
        void copyData();

        vsg::ref_ptr<CloudsBinding> loadData();

    private:
        vsg::ref_ptr<vsg::Device> _device;
        vsg::ref_ptr<vsg::PhysicalDevice> _physicalDevice;

        vsg::ref_ptr<AtmosphereModelSettings> _settings;

        vsg::ref_ptr<vsg::ShaderStage> _shapeShader;
        vsg::ref_ptr<vsg::ShaderStage> _detailShader;

        vsg::ref_ptr<Image> _shapeNoiseTexture;
        vsg::ref_ptr<Image> _detailNoiseTexture;
        vsg::ref_ptr<Image> _curlNoiseTexture;
        vsg::ref_ptr<Image> _blueNoiseTexture;
    };
    extern vsg::ref_ptr<CloudsBinding> createCloudsData(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Options> options, vsg::ref_ptr<AtmosphereModelSettings> settings);
}
#endif // CLOUDS_H
