#include "BRDF.h"
#include <vsg/commands/Commands.h>
#include <vsg/commands/Dispatch.h>
#include <vsg/state/ComputePipeline.h>
#include <vsg/app/CompileTraversal.h>
#include <vsg/vk/SubmitCommands.h>

namespace atmosphere {

    BRDFGenerator::BRDFGenerator(vsg::ref_ptr<vsg::Device> device,
                                     vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice,
                                     vsg::ref_ptr<AtmosphereModelSettings> settings,
                                     vsg::ref_ptr<vsg::Options> options)
        : _device(device)
        , _physicalDevice(physicalDevice)
        , _settings(settings)
    {
        auto sampler = vsg::Sampler::create();
        _BRDFlutTexture = Image::create(VkExtent3D{settings->BRDFlutSize, settings->BRDFlutSize, 1}, sampler, VK_FORMAT_R16G16_SFLOAT);
        _BRDFlutTexture->allocateTexture(_device);

        _lutShader = vsg::ShaderStage::create(VK_SHADER_STAGE_COMPUTE_BIT, "main", options->getRefObject<vsg::ShaderModule>(std::string(BRDFlutShader)));

        _lutShader->module->hints = vsg::ShaderCompileSettings::create();
        _lutShader->module->hints->generateDebugInfo = true;
    }

    void BRDFGenerator::initialize()
    {
        auto noise = vsg::Commands::create();

        // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}
        };
        auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{descriptorSetLayout}, vsg::PushConstantRanges{});

        auto buffer = vsg::DescriptorImage::create(_BRDFlutTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

        _lutShader->specializationConstants.insert_or_assign(1, vsg::intValue::create(_settings->numThreads));
        _lutShader->specializationConstants.insert_or_assign(2, vsg::intValue::create(_settings->BRDFsamples));
        _lutShader->specializationConstants.insert_or_assign(3, vsg::intValue::create(_settings->BRDFlutSize));

        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, vsg::Descriptors{buffer});
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, descriptorSet);

        auto imageBarrier = _BRDFlutTexture->transitionWriteBarrier();
        auto transitionBarrier = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, imageBarrier);

        auto pipeline = vsg::ComputePipeline::create(pipelineLayout, _lutShader);
        auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
        noise->addChild(transitionBarrier);
        noise->addChild(bindPipeline);
        noise->addChild(bindTextures);
        auto size = uint32_t(ceil(float(_settings->BRDFlutSize) / float(_settings->numThreads)));
        noise->addChild(vsg::Dispatch::create(size, size, 1));

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

    void BRDFGenerator::copyData()
    {
        auto copy = _BRDFlutTexture->copyData(_device);

        auto fence = vsg::Fence::create(_device);
        auto queueFamilyIndex = _physicalDevice->getQueueFamily(VK_QUEUE_GRAPHICS_BIT);
        auto commandPool = vsg::CommandPool::create(_device, queueFamilyIndex);
        auto queue = _device->getQueue(queueFamilyIndex);

        vsg::submitCommandsToQueue(commandPool, fence, 100000000000, queue, [&](vsg::CommandBuffer& commandBuffer) {
            copy->record(commandBuffer);
        });

        _BRDFlutTexture->mapData(_device);
    }

    vsg::ref_ptr<BRDFBinding> BRDFGenerator::loadData()
    {
        auto binding = BRDFBinding::create(_settings->prefilteredSize, _settings->prefilteredLayers);
        binding->BRDFlutTexture = _BRDFlutTexture;

        return binding;
    }
}
