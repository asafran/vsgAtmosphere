#include "AtmosphereImage.h"
#include <vsg/commands/Commands.h>
#include <vsg/state/ComputePipeline.h>
#include <vsg/state/DescriptorSetLayout.h>
#include <vsg/state/PipelineLayout.h>
#include <vsg/vk/Fence.h>

namespace atmosphere {

    Image::Image(VkExtent3D extent, vsg::ref_ptr<vsg::Sampler> sampler, VkFormat format)
    {
        vsg::ref_ptr<vsg::Image> image = vsg::Image::create();
        image->usage |= (VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT);
        image->format = format;
        image->mipLevels = 1;
        image->extent = extent;
        image->imageType = extent.depth > 1 ? VK_IMAGE_TYPE_3D : VK_IMAGE_TYPE_2D;
        image->arrayLayers = 1;

        auto imageView = vsg::ImageView::create(image, VK_IMAGE_ASPECT_COLOR_BIT);
        if(!sampler)
            sampler = vsg::Sampler::create();
        imageInfo = vsg::ImageInfo::create(sampler, imageView, VK_IMAGE_LAYOUT_GENERAL);
    }

    Image::Image(vsg::ref_ptr<vsg::Data> in_data, vsg::ref_ptr<vsg::Sampler> sampler)
        : imageInfo(vsg::ImageInfo::create(sampler, in_data))
        , data(in_data)
    {

    }

    Image::Image()
    {

    }

    Image::~Image()
    {

    }

    int Image::compare(const Object &rhs_object) const
    {
        int result = Object::compare(rhs_object);
        if (result != 0) return result;

        auto& rhs = static_cast<decltype(*this)>(rhs_object);
        return compare_pointer(data, rhs.data);
    }

    void Image::read(vsg::Input &input)
    {
        vsg::ref_ptr<vsg::Sampler> sampler;
        input.readObject("sampler", sampler);
        input.readObject("data", data);

        imageInfo = vsg::ImageInfo::create(sampler, data);
    }

    void Image::write(vsg::Output &output) const
    {
        output.writeObject("sampler", imageInfo->sampler);
        output.writeObject("data", data);
    }

    vsg::ref_ptr<vsg::ImageMemoryBarrier> Image::transitionWriteBarrier()
    {
        if(!_transitionWriteBarrier)
        {
            _transitionWriteBarrier = vsg::ImageMemoryBarrier::create();
            _transitionWriteBarrier->srcAccessMask = 0;
            _transitionWriteBarrier->dstAccessMask = VK_ACCESS_MEMORY_WRITE_BIT;
            _transitionWriteBarrier->oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
            _transitionWriteBarrier->newLayout = VK_IMAGE_LAYOUT_GENERAL;
            _transitionWriteBarrier->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionWriteBarrier->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionWriteBarrier->image = imageInfo->imageView->image;
            _transitionWriteBarrier->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            _transitionWriteBarrier->subresourceRange.baseArrayLayer = 0;
            _transitionWriteBarrier->subresourceRange.layerCount = 1;
            _transitionWriteBarrier->subresourceRange.levelCount = 1;
            _transitionWriteBarrier->subresourceRange.baseMipLevel = 0;
        }
        return _transitionWriteBarrier;
    }

    vsg::ref_ptr<vsg::ImageMemoryBarrier> Image::transitionReadBarrier()
    {
        if(!_transitionReadBarrier)
        {
            _transitionReadBarrier = vsg::ImageMemoryBarrier::create();
            _transitionReadBarrier->srcAccessMask = 0;
            _transitionReadBarrier->dstAccessMask = VK_ACCESS_MEMORY_READ_BIT;
            _transitionReadBarrier->oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
            _transitionReadBarrier->newLayout = VK_IMAGE_LAYOUT_GENERAL;
            _transitionReadBarrier->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionReadBarrier->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionReadBarrier->image = imageInfo->imageView->image;
            _transitionReadBarrier->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            _transitionReadBarrier->subresourceRange.baseArrayLayer = 0;
            _transitionReadBarrier->subresourceRange.layerCount = 1;
            _transitionReadBarrier->subresourceRange.levelCount = 1;
            _transitionReadBarrier->subresourceRange.baseMipLevel = 0;
        }
        return _transitionReadBarrier;
    }

    vsg::ref_ptr<vsg::ImageMemoryBarrier> Image::transitionOptimalBarrier()
    {
        if(!_transitionOptimal)
        {
            _transitionOptimal = vsg::ImageMemoryBarrier::create();
            _transitionOptimal->srcAccessMask = 0;
            _transitionOptimal->dstAccessMask = VK_ACCESS_MEMORY_READ_BIT;
            _transitionOptimal->oldLayout = VK_IMAGE_LAYOUT_GENERAL;
            _transitionOptimal->newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
            _transitionOptimal->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionOptimal->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _transitionOptimal->image = imageInfo->imageView->image;
            _transitionOptimal->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            _transitionOptimal->subresourceRange.baseArrayLayer = 0;
            _transitionOptimal->subresourceRange.layerCount = 1;
            _transitionOptimal->subresourceRange.levelCount = 1;
            _transitionOptimal->subresourceRange.baseMipLevel = 0;
        }
        return _transitionOptimal;
    }

    vsg::ref_ptr<vsg::ImageMemoryBarrier> Image::writeBarrier()
    {
        if(!_writeBarrier)
        {
            _writeBarrier = vsg::ImageMemoryBarrier::create();
            _writeBarrier->srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT ;
            _writeBarrier->dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
            _writeBarrier->oldLayout = VK_IMAGE_LAYOUT_GENERAL;
            _writeBarrier->newLayout = VK_IMAGE_LAYOUT_GENERAL;
            _writeBarrier->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _writeBarrier->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
            _writeBarrier->image = imageInfo->imageView->image;
            _writeBarrier->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
            _writeBarrier->subresourceRange.baseArrayLayer = 0;
            _writeBarrier->subresourceRange.layerCount = 1;
            _writeBarrier->subresourceRange.levelCount = 1;
            _writeBarrier->subresourceRange.baseMipLevel = 0;
        }
        return _writeBarrier;
    }

    void Image::allocateTexture(vsg::Device* device, bool init)
    {
        auto &image = imageInfo->imageView->image;
        if(init)
            switch (image->imageType) {
            case VK_IMAGE_TYPE_2D:
                image->data = vsg::vec4Array2D::create(image->extent.width, image->extent.height, vsg::vec4{0.0f, 0.0f, 0.0f, 1.0f}, vsg::Data::Properties{image->format});
                break;
            case VK_IMAGE_TYPE_3D:
                image->data = vsg::vec4Array3D::create(image->extent.width, image->extent.height, image->extent.depth, vsg::vec4{0.0f, 0.0f, 0.0f, 1.0f}, vsg::Data::Properties{image->format});
            default:
                break;
            }


        image->compile(device);
        image->allocateAndBindMemory(device, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    }

    vsg::ref_ptr<vsg::CopyImage> Image::copyData(vsg::Device *device)
    {
        auto sourceImage = imageInfo->imageView->image;

        auto width = sourceImage->extent.width;
        auto height = sourceImage->extent.height;
        auto depth = sourceImage->extent.depth;

        VkFormat sourceImageFormat = sourceImage->format;
        VkFormat targetImageFormat = sourceImageFormat;

        _converted = vsg::Image::create();
        _converted->imageType = sourceImage->imageType;
        _converted->format = targetImageFormat;
        _converted->extent.width = width;
        _converted->extent.height = height;
        _converted->extent.depth = depth;
        _converted->arrayLayers = 1;
        _converted->mipLevels = 1;
        _converted->initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        _converted->samples = VK_SAMPLE_COUNT_1_BIT;
        _converted->tiling = VK_IMAGE_TILING_LINEAR;
        _converted->usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT;

        _converted->compile(device);

        auto deviceMemory = vsg::DeviceMemory::create(device, _converted->getMemoryRequirements(device->deviceID), VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);

        _converted->bind(deviceMemory, 0);

        VkImageCopy region{};
        region.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        region.srcSubresource.layerCount = 1;
        region.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        region.dstSubresource.layerCount = 1;
        region.extent.width = width;
        region.extent.height = height;
        region.extent.depth = depth;

        auto copyImage = vsg::CopyImage::create();
        copyImage->srcImage = sourceImage;
        copyImage->srcImageLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
        copyImage->dstImage = _converted;
        copyImage->dstImageLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        copyImage->regions.push_back(region);

        return copyImage;
    }

    void Image::mapData(vsg::Device *device)
    {
        VkImageSubresource subResource{VK_IMAGE_ASPECT_COLOR_BIT, 0, 0};
        VkSubresourceLayout subResourceLayout;
        vkGetImageSubresourceLayout(*device, _converted->vk(device->deviceID), &subResource, &subResourceLayout);

        auto deviceMemory = _converted->getDeviceMemory(device->deviceID);
        auto sourceImage = imageInfo->imageView->image;
        auto width = sourceImage->extent.width;
        auto height = sourceImage->extent.height;
        auto depth = sourceImage->extent.depth;

        VkFormat targetImageFormat =  _converted->format;

        size_t destRowWidth = width * sizeof(vsg::vec4);
        if (destRowWidth == subResourceLayout.rowPitch)
        {
            if(sourceImage->imageType == VK_IMAGE_TYPE_3D)
                data = vsg::MappedData<vsg::vec4Array3D>::create(deviceMemory, subResourceLayout.offset, 0, vsg::Data::Properties{targetImageFormat}, width, height, depth); // deviceMemory, offset, flags and dimensions
            else
                data = vsg::MappedData<vsg::vec4Array2D>::create(deviceMemory, subResourceLayout.offset, 0, vsg::Data::Properties{targetImageFormat}, width, height); // deviceMemory, offset, flags and dimensions
        }
        else
        {
            if(sourceImage->imageType != VK_IMAGE_TYPE_3D)
            {
                // Map the buffer memory and assign as a ubyteArray that will automatically unmap itself on destruction.
                // A ubyteArray is used as the graphics buffer memory is not contiguous like vsg::Array2D, so map to a flat buffer first then copy to Array2D.
                auto mappedData = vsg::MappedData<vsg::floatArray>::create(deviceMemory, subResourceLayout.offset, 0, vsg::Data::Properties{targetImageFormat}, subResourceLayout.rowPitch*height);
                data = vsg::vec4Array2D::create(width, height, vsg::Data::Properties{targetImageFormat});
                for (uint32_t row = 0; row < height; ++row)
                {
                    std::memcpy(data->dataPointer(row*width), mappedData->dataPointer(row * subResourceLayout.rowPitch), destRowWidth);
                }
            }
        }
        _converted = {};
    }

    /*vsg::ref_ptr<vsg::Commands> AtmosphereImage::generateTextureCommands(const vsg::DescriptorSetLayouts &dsls, int numThreads)
    {
        auto commands = vsg::Commands::create();
        auto pipelineLayout = vsg::PipelineLayout::create(dsls, vsg::PushConstantRanges{});
        auto pipeline = vsg::ComputePipeline::create(pipelineLayout, genShaderStage);
        auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
        commands->addChild(bindPipeline);

        auto &extent = imageInfo->imageView->image->extent;
        auto width = uint32_t(ceil(float(extent.width) / float(numThreads)));
        auto height = uint32_t(ceil(float(extent.height) / float(numThreads)));
        auto depth = uint32_t(ceil(float(extent.depth) / float(numThreads)));
        commands->addChild(vsg::Dispatch::create(width, height, depth));

        return commands;
    }*/
}
