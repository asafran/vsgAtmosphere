#include "Clouds.h"
#include "AtmoshpereConstants.h"
#include <vsg/commands/Commands.h>
#include <vsg/commands/Dispatch.h>
#include <vsg/state/ComputePipeline.h>
#include <vsg/app/CompileTraversal.h>
#include <vsg/vk/SubmitCommands.h>

namespace atmosphere {

CloudsGenerator::CloudsGenerator(vsg::ref_ptr<vsg::Device> device,
                                 vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice,
                                 vsg::ref_ptr<AtmosphereModelSettings> settings,
                                 vsg::ref_ptr<vsg::Options> options)
        : _device(device)
        , _physicalDevice(physicalDevice)
        , _settings(settings)
    {
        auto sampler = vsg::Sampler::create();
        _shapeNoiseTexture = Image::create(VkExtent3D{_settings->shapeNoiseSize, _settings->shapeNoiseSize, _settings->shapeNoiseSize}, sampler, VK_FORMAT_R16G16B16A16_SFLOAT);
        _shapeNoiseTexture->allocateTexture(_device);

        _detailNoiseTexture = Image::create(VkExtent3D{_settings->detailNoiseSize, _settings->detailNoiseSize, _settings->detailNoiseSize}, sampler, VK_FORMAT_R16G16B16A16_SFLOAT);
        _detailNoiseTexture->allocateTexture(_device);

        _curlNoiseTexture = Image::create(options->getRefObject<vsg::Data>(std::string(curlNoiseTextureName)), sampler);
        _blueNoiseTexture = Image::create(options->getRefObject<vsg::Data>(std::string(blueNoiseTextureName)), sampler);

        _detailShader = vsg::ShaderStage::create(VK_SHADER_STAGE_COMPUTE_BIT, "main", options->getRefObject<vsg::ShaderModule>(std::string(detailNoiseShader)));
        _shapeShader = vsg::ShaderStage::create(VK_SHADER_STAGE_COMPUTE_BIT, "main", options->getRefObject<vsg::ShaderModule>(std::string(shapeNoiseShader)));
    }

    void CloudsGenerator::initialize()
    {
        auto noise = vsg::Commands::create();

        // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
            {1, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}
        };
        auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{descriptorSetLayout}, vsg::PushConstantRanges{});

        auto shapeBuffer = vsg::DescriptorImage::create(_shapeNoiseTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
        auto detailBuffer = vsg::DescriptorImage::create(_detailNoiseTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

        _detailShader->specializationConstants.insert_or_assign(0, vsg::intValue::create(_settings->numThreads));
        _detailShader->specializationConstants.insert_or_assign(1, vsg::intValue::create(_settings->detailNoiseSize));
        _shapeShader->specializationConstants.insert_or_assign(0, vsg::intValue::create(_settings->numThreads));
        _shapeShader->specializationConstants.insert_or_assign(1, vsg::intValue::create(_settings->shapeNoiseSize));

        vsg::Descriptors descriptors{shapeBuffer, detailBuffer};
        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, descriptorSet);
        {

            auto pipeline = vsg::ComputePipeline::create(pipelineLayout, _shapeShader);
            auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
            noise->addChild(bindPipeline);
            noise->addChild(bindTextures);
            auto size = uint32_t(ceil(float(_settings->shapeNoiseSize) / float(_settings->numThreads)));
            noise->addChild(vsg::Dispatch::create(size, size, size));
        }

        {
            auto pipeline = vsg::ComputePipeline::create(pipelineLayout, _detailShader);
            auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
            noise->addChild(bindPipeline);
            noise->addChild(bindTextures);
            auto size = uint32_t(ceil(float(_settings->detailNoiseSize) / float(_settings->numThreads)));
            noise->addChild(vsg::Dispatch::create(size, size, size));
        }

        // compile the Vulkan objects
        auto compileTraversal = vsg::CompileTraversal::create(_device);
        auto context = compileTraversal->contexts.front();
        noise->accept(*compileTraversal);

        int computeQueueFamily = _physicalDevice->getQueueFamily(VK_QUEUE_COMPUTE_BIT);

        context->commandPool = vsg::CommandPool::create(_device, computeQueueFamily);

        auto fence = vsg::Fence::create(_device);
        auto computeQueue = _device->getQueue(computeQueueFamily);

        vsg::submitCommandsToQueue(context->commandPool, fence, 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
            noise->record(commandBuffer);
        });
    }

    void CloudsGenerator::copyData()
    {
        auto commands = vsg::Commands::create();
        commands->addChild(_shapeNoiseTexture->copyData(_device));
        commands->addChild(_detailNoiseTexture->copyData(_device));

        auto fence = vsg::Fence::create(_device);
        auto queueFamilyIndex = _physicalDevice->getQueueFamily(VK_QUEUE_GRAPHICS_BIT);
        auto commandPool = vsg::CommandPool::create(_device, queueFamilyIndex);
        auto queue = _device->getQueue(queueFamilyIndex);

        vsg::submitCommandsToQueue(commandPool, fence, 100000000000, queue, [&](vsg::CommandBuffer& commandBuffer) {
            commands->record(commandBuffer);
        });

        _shapeNoiseTexture->mapData(_device);
        _detailNoiseTexture->mapData(_device);

    }

    vsg::ref_ptr<CloudsBinding> CloudsGenerator::loadData()
    {
        auto cloudsBinding = CloudsBinding::create();
        cloudsBinding->curlNoiseTexture = _curlNoiseTexture;
        cloudsBinding->blueNoiseTexture = _blueNoiseTexture;
        cloudsBinding->shapeNoiseTexture = _shapeNoiseTexture;
        cloudsBinding->detailNoiseTexture = _detailNoiseTexture;

        return cloudsBinding;
    }
}
