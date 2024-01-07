#ifndef TASK_H
#define TASK_H

#include <vsg/state/ComputePipeline.h>
#include "AtmosphereBinding.h"

namespace atmosphere {

    class Task : public vsg::Inherit<vsg::Object, Task>
    {
    public:

        vsg::ref_ptr<vsg::Commands> createTaskCommands();
        vsg::ref_ptr<vsg::DescriptorSetLayout> createTexturesSetLayout();
        vsg::ref_ptr<vsg::DescriptorSet> createTexturesSet(vsg::ref_ptr<vsg::DescriptorSetLayout> descriptorSetLayout);

        std::vector<vsg::ref_ptr<Image>> readImages;
        std::vector<vsg::ref_ptr<Image>> writeImages;

        int32_t numThreads = 32;

        vsg::ref_ptr<ComputeParametersBinding> parameters;
        vsg::ref_ptr<vsg::ShaderStage> shader;

    private:
        vsg::ref_ptr<vsg::Commands> _taskCommands;
    };
}



#endif // TASK_H
