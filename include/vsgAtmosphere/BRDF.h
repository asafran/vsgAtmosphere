#ifndef BRDF_H
#define BRDF_H

#include "Atmosphere.h"
#include "AtmosphereBinding.h"
#include "AtmosphereImage.h"
#include <vsg/app/Window.h>
#include <vsg/vk/Device.h>
#include <vsg/state/ImageInfo.h>
#include <vsg/io/Options.h>

namespace atmosphere {

    class BRDFGenerator : public vsg::Inherit<vsg::Object, BRDFGenerator>
    {
    public:
        BRDFGenerator(vsg::ref_ptr<vsg::Device> device,
                        vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice,
                        vsg::ref_ptr<AtmosphereModelSettings> settings,
                        vsg::ref_ptr<vsg::Options> options);

        void initialize();
        void copyData();

        vsg::ref_ptr<BRDFBinding> loadData();

    private:
        vsg::ref_ptr<vsg::Device> _device;
        vsg::ref_ptr<vsg::PhysicalDevice> _physicalDevice;

        vsg::ref_ptr<AtmosphereModelSettings> _settings;

        vsg::ref_ptr<vsg::ShaderStage> _lutShader;

        vsg::ref_ptr<Image> _BRDFlutTexture;
    };
}

#endif // BRDF_H
