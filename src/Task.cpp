#include "Task.h"
#include <vsg/state/BindDescriptorSet.h>
#include <vsg/state/DescriptorImage.h>
#include <vsg/commands/Commands.h>
#include <vsg/commands/Dispatch.h>
#include <vsg/io/Options.h>
#include <vsg/io/Logger.h>

namespace atmosphere {

    vsg::ref_ptr<vsg::Commands> Task::createTaskCommands()
    {
        if(!_taskCommands)
        {
            _taskCommands = vsg::Commands::create();

            vsg::PushConstantRanges orderPc{
                {VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(int)}
            };
            auto parametersLayout = parameters->descriptorSetLayout;
            auto texturesLayout = createTexturesSetLayout();
            auto texturesSet = createTexturesSet(texturesLayout);
            auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesLayout, parametersLayout}, orderPc);
            auto pipeline = vsg::ComputePipeline::create(pipelineLayout, shader);
            auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
            auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, TEXTURES_DESCRIPTOR_SET, texturesSet);
            auto bindParameters = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, PARAMETERS_DESCRIPTOR_SET, parameters->descriptorSet);
            _taskCommands->addChild(bindPipeline);
            _taskCommands->addChild(bindTextures);
            _taskCommands->addChild(bindParameters);

            auto extent = writeImages.front()->imageInfo->imageView->image->extent;
            _taskCommands->addChild(vsg::Dispatch::create(uint32_t(ceil(float(extent.width) / float(numThreads))),
                                                          uint32_t(ceil(float(extent.height) / float(numThreads))),
                                                          uint32_t(ceil(float(extent.depth) / float(numThreads)))));
        }
        return _taskCommands;
    }

    vsg::ref_ptr<vsg::DescriptorSetLayout> Task::createTexturesSetLayout()
    {
        vsg::DescriptorSetLayoutBindings descriptorBindings;
        uint32_t binding = 0;
        for (auto image : writeImages)
        {
            descriptorBindings.push_back({binding, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr});
            binding++;
        }
        for (auto image : readImages)
        {
            descriptorBindings.push_back({binding, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr});
            binding++;
        }
        return vsg::DescriptorSetLayout::create(descriptorBindings);
    }

    vsg::ref_ptr<vsg::DescriptorSet> Task::createTexturesSet(vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout)
    {
        vsg::Descriptors descriptors;
        uint32_t binding = 0;
        for (auto image : writeImages)
        {
            descriptors.push_back(vsg::DescriptorImage::create(image->imageInfo, binding, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE));
            binding++;
        }
        for (auto image : readImages)
        {
            descriptors.push_back(vsg::DescriptorImage::create(image->imageInfo, binding, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER));
            binding++;
        }
        return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
    }
}
