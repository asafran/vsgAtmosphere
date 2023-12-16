#ifndef ATMOSPHEREIMAGE_H
#define ATMOSPHEREIMAGE_H

#include <vsg/commands/CopyImage.h>
#include <vsg/state/DescriptorSetLayout.h>
#include <vsg/state/ImageInfo.h>

namespace atmosphere {

class Image : public vsg::Inherit<vsg::Object, Image>
    {
    public:
        Image(VkExtent3D extent, vsg::ref_ptr<vsg::Sampler> sampler = {}, VkFormat format = VK_FORMAT_R32G32B32A32_SFLOAT);
<<<<<<< Updated upstream
        Image();
        virtual ~Image();

=======
        Image(vsg::ref_ptr<vsg::Data> in_data, vsg::ref_ptr<vsg::Sampler> sampler = {});
        Image();
        virtual ~Image();

        int compare(const Object& rhs_object) const override;

        void read(vsg::Input& input) override;
        void write(vsg::Output& output) const override;

>>>>>>> Stashed changes
        vsg::ref_ptr<vsg::ImageInfo> imageInfo;
        vsg::ref_ptr<vsg::Data> data;

        void allocateTexture(vsg::Device *device, bool init = false);
        vsg::ref_ptr<vsg::CopyImage> copyData(vsg::Device *device);
        void mapData(vsg::Device *device);

    private:
        vsg::ref_ptr<vsg::Image> _converted;
    };

}
<<<<<<< Updated upstream
=======

EVSG_type_name(atmosphere::Image)

>>>>>>> Stashed changes
#endif // ATMOSPHEREIMAGE_H
