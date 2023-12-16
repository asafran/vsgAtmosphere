#include "Atmosphere.h"
#include "AtmosphereLighting.h"
#include "AtmosphereTools.h"
#include "Clouds.h"
#include "InverseMatrices.h"

#include <vsg/all.h>


namespace atmosphere {

void AtmosphereGenerator::copyData()
{
    auto commands = vsg::Commands::create();
    commands->addChild(_transmittanceTexture->copyData(_device));
    commands->addChild(_irradianceTexture->copyData(_device));
    commands->addChild(_scatteringTexture->copyData(_device));
    commands->addChild(_singleMieScatteringTexture->copyData(_device));

    auto fence = vsg::Fence::create(_device);
    auto queueFamilyIndex = _physicalDevice->getQueueFamily(VK_QUEUE_GRAPHICS_BIT);
    auto commandPool = vsg::CommandPool::create(_device, queueFamilyIndex);
    auto queue = _device->getQueue(queueFamilyIndex);

    vsg::submitCommandsToQueue(commandPool, fence, 100000000000, queue, [&](vsg::CommandBuffer& commandBuffer) {
        commands->record(commandBuffer);
    });

    _transmittanceTexture->mapData(_device);
    _irradianceTexture->mapData(_device);
    _scatteringTexture->mapData(_device);
    _singleMieScatteringTexture->mapData(_device);
}

/*
vsg::ref_ptr<Clouds> loadClouds(const vsg::Path &path, vsg::ref_ptr<const vsg::Options> options)
{

}
*/
AtmosphereModelSettings::AtmosphereModelSettings(vsg::ref_ptr<vsg::EllipsoidModel> model)
    : ellipsoidModel(model)
{
}

AtmosphereModelSettings::AtmosphereModelSettings()
{

}

AtmosphereModelSettings::~AtmosphereModelSettings()
{

}

void AtmosphereModelSettings::read(vsg::Input &input)
{
    input.read("numThreads", numThreads);

    input.read("transmittanceWidth", transmittanceWidth);
    input.read("transmittanceHeight", transmittanceHeight);

    input.read("scaterringR", scaterringR);
    input.read("scaterringMU", scaterringMU);
    input.read("scaterringMU_S", scaterringMU_S);
    input.read("scaterringNU", scaterringNU);

    input.read("irradianceWidth", irradianceWidth);
    input.read("irradianceHeight", irradianceHeight);

    input.read("precomputedWavelenghts", precomputedWavelenghts);

    input.read("transmittanceHeight", transmittanceHeight);

    input.read("sunAngularRadius", sunAngularRadius);

    input.read("rayleighDensityLayer", rayleighDensityLayer);
    input.read("mieDensityLayer", mieDensityLayer);

    input.read("absorptionDensityLayer0", absorptionDensityLayer0);
    input.read("absorptionDensityLayer1", absorptionDensityLayer1);

    input.read("miePhaseFunction_g", miePhaseFunction_g);

    input.read("maxSunZenithAngle", maxSunZenithAngle);
    input.read("lengthUnitInMeters", lengthUnitInMeters);

    input.read("ellipsoidModel", ellipsoidModel);

    input.read("kSolarIrradiance", kSolarIrradiance);
    input.read("kOzoneCrossSection", kOzoneCrossSection);

    input.read("kDobsonUnit", kDobsonUnit);
    input.read("kMaxOzoneNumberDensity", kMaxOzoneNumberDensity);
    input.read("kConstantSolarIrradiance", kConstantSolarIrradiance);

    input.read("kAtmoshpereHeight", atmoshpereHeight);

    input.read("kRayleigh", kRayleigh);
    input.read("kRayleighScaleHeight", kRayleighScaleHeight);
    input.read("kMieScaleHeight", kMieScaleHeight);
    input.read("kMieAngstromAlpha", kMieAngstromAlpha);
    input.read("kMieAngstromBeta", kMieAngstromBeta);
    input.read("kMieSingleScatteringAlbedo", kMieSingleScatteringAlbedo);
    input.read("kGroundAlbedo", kGroundAlbedo);

    input.read("use_constant_solar_spectrum_", transmittanceHeight);
    input.read("use_ozone", use_ozone);
}

void AtmosphereModelSettings::write(vsg::Output &output) const
{
    output.write("numThreads", numThreads);

    output.write("transmittanceWidth", transmittanceWidth);
    output.write("transmittanceHeight", transmittanceHeight);

    output.write("scaterringR", scaterringR);
    output.write("scaterringMU", scaterringMU);
    output.write("scaterringMU_S", scaterringMU_S);
    output.write("scaterringNU", scaterringNU);

    output.write("irradianceWidth", irradianceWidth);
    output.write("irradianceHeight", irradianceHeight);

    output.write("precomputedWavelenghts", precomputedWavelenghts);

    output.write("transmittanceHeight", transmittanceHeight);

    output.write("sunAngularRadius", sunAngularRadius);

    output.write("rayleighDensityLayer", rayleighDensityLayer);
    output.write("mieDensityLayer", mieDensityLayer);

    output.write("absorptionDensityLayer0", absorptionDensityLayer0);
    output.write("absorptionDensityLayer1", absorptionDensityLayer1);

    output.write("miePhaseFunction_g", miePhaseFunction_g);

    output.write("maxSunZenithAngle", maxSunZenithAngle);
    output.write("lengthUnitInMeters", lengthUnitInMeters);

    output.write("ellipsoidModel", ellipsoidModel);

    output.write("kSolarIrradiance", kSolarIrradiance);
    output.write("kOzoneCrossSection", kOzoneCrossSection);

    output.write("kDobsonUnit", kDobsonUnit);
    output.write("kMaxOzoneNumberDensity", kMaxOzoneNumberDensity);
    output.write("kConstantSolarIrradiance", kConstantSolarIrradiance);

    output.write("kAtmoshpereHeight", atmoshpereHeight);

    output.write("kRayleigh", kRayleigh);
    output.write("kRayleighScaleHeight", kRayleighScaleHeight);
    output.write("kMieScaleHeight", kMieScaleHeight);
    output.write("kMieAngstromAlpha", kMieAngstromAlpha);
    output.write("kMieAngstromBeta", kMieAngstromBeta);
    output.write("kMieSingleScatteringAlbedo", kMieSingleScatteringAlbedo);
    output.write("kGroundAlbedo", kGroundAlbedo);

    output.write("use_constant_solar_spectrum_", transmittanceHeight);
    output.write("use_ozone", use_ozone);
}

vsg::ShaderStage::SpecializationConstants AtmosphereModelSettings::getComputeConstants() const
{
    return vsg::ShaderStage::SpecializationConstants{
        {0, vsg::intValue::create(numThreads)},
        //{1, vsg::intValue::create(numViewerThreads)},
        //{2, vsg::intValue::create(cubeSize)},

        {10, vsg::intValue::create(transmittanceWidth)},
        {11, vsg::intValue::create(transmittanceHeight)},

        {12, vsg::intValue::create(scaterringR)},
        {13, vsg::intValue::create(scaterringMU)},
        {14, vsg::intValue::create(scaterringMU_S)},
        {15, vsg::intValue::create(scaterringNU)},

        {16, vsg::intValue::create(irradianceWidth)},
        {17, vsg::intValue::create(irradianceHeight)},

        {36, vsg::floatValue::create(sunAngularRadius)},
        {37, vsg::floatValue::create(static_cast<float>(ellipsoidModel->radiusEquator() / lengthUnitInMeters))},
        {38, vsg::floatValue::create(static_cast<float>((ellipsoidModel->radiusEquator() + atmoshpereHeight) / lengthUnitInMeters))},
        {39, vsg::floatValue::create(miePhaseFunction_g)},
        {40, vsg::floatValue::create(static_cast<float>(std::cos(maxSunZenithAngle)))}
    };
}

// -----------------------------------------------------------------------------------------------------------------------------------

AtmosphereGenerator::AtmosphereGenerator(vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Device> device, vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice, vsg::ref_ptr<vsg::Options> options)
    : _settings(settings)
    , _device(device)
    , _physicalDevice(physicalDevice)
    , _options(options)
{
    compileSettings = vsg::ShaderCompileSettings::create();

    constexpr int kLambdaMinN = static_cast<int>(kLambdaMin);
    constexpr int kLambdaMaxN = static_cast<int>(kLambdaMax);
/*
    numThreads = settings->numThreads;

    transmittanceWidth = settings->transmittanceWidth;
    transmittanceHeight = settings->transmittanceHeight;

    scaterringR = settings->scaterringR;
    scaterringMU = settings->scaterringMU;
    scaterringMU_S = settings->scaterringMU_S;
    scaterringNU = settings->scaterringNU;

    irradianceWidth = settings->irradianceWidth;
    irradianceHeight = settings->irradianceHeight;
*/


    for (int l = kLambdaMinN; l <= kLambdaMaxN; l += 10) {
        double lambda = static_cast<double>(l) * 1e-3;  // micro-meters
        double mie =
                settings->kMieAngstromBeta / settings->kMieScaleHeight * pow(lambda, -settings->kMieAngstromAlpha);
        waveLengths.push_back(l);
        if (settings->use_constant_solar_spectrum_) {
            solarIrradiance.push_back(settings->kConstantSolarIrradiance);
        } else {
            solarIrradiance.push_back(settings->kSolarIrradiance[(l - kLambdaMinN) / 10]);
        }
        rayleighScattering.push_back(settings->kRayleigh * pow(lambda, -4));
        mieScattering.push_back(mie * settings->kMieSingleScatteringAlbedo);
        mieExtinction.push_back(mie);
        absorptionExtinction.push_back(settings->use_ozone ?
                                            settings->kMaxOzoneNumberDensity * settings->kOzoneCrossSection[(l - kLambdaMinN) / 10] :
                                        0.0);
        groundAlbedo.push_back(settings->kGroundAlbedo);
    }
}
AtmosphereGenerator::AtmosphereGenerator(vsg::ref_ptr<vsg::Device> device, vsg::ref_ptr<vsg::PhysicalDevice> physicalDevice, vsg::ref_ptr<vsg::Options> options)
    : _device(device)
    , _physicalDevice(physicalDevice)
    , _options(options)
{
    compileSettings = vsg::ShaderCompileSettings::create();
}

// -----------------------------------------------------------------------------------------------------------------------------------

AtmosphereGenerator::~AtmosphereGenerator()
{
}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereGenerator::initialize()
{
    generateTextures();

    _computeConstants = _settings->getComputeConstants();// assignComputeConstants();
    assignRenderConstants();

    auto memoryBarrier = vsg::MemoryBarrier::create();
    memoryBarrier->srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    memoryBarrier->dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    auto memoryBarrierCmd = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, memoryBarrier);

    auto singlePass = vsg::Commands::create();
    auto transmittanceCommands = vsg::Commands::create();
    auto multipleScattering = vsg::Commands::create();
    auto lfr = vsg::mat4Value::create(vsg::mat4());
    auto parameters = vsg::Value<Parameters>::create(Parameters());
    auto profiles = vsg::Array<DensityProfileLayer>::create(4);
    profiles->set(0, _settings->rayleighDensityLayer);
    profiles->set(1, _settings->mieDensityLayer);
    profiles->set(2, _settings->absorptionDensityLayer0);
    profiles->set(3, _settings->absorptionDensityLayer1);

    auto pLayout = parametersLayout();

    auto lfrBuffer = vsg::DescriptorBuffer::create(lfr, 0, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER);
    auto parametersBuffer = vsg::DescriptorBuffer::create(parameters, 1, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER);
    auto proflesArray = vsg::DescriptorBuffer::create(profiles, 2, 0, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER);

    vsg::ref_ptr<vsg::BindDescriptorSet> bindParameters;
    vsg::ref_ptr<vsg::BindDescriptorSet> bindOrder;

    {
        vsg::Descriptors descriptors{lfrBuffer, parametersBuffer, proflesArray};
        auto descriptorSet = vsg::DescriptorSet::create(pLayout, descriptors);
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{pLayout}, vsg::PushConstantRanges{});
        bindParameters = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, 1, descriptorSet);
        bindParameters->slot = 2;
    }

    auto oLayout = orderLayout();

    auto order = vsg::intValue::create(0);
    auto orderBuffer = vsg::DescriptorBuffer::create(order, 0, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER);

    {
        vsg::Descriptors descriptors{orderBuffer};
        auto descriptorSet = vsg::DescriptorSet::create(oLayout, descriptors);
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{oLayout}, vsg::PushConstantRanges{});
        bindOrder = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, 2, descriptorSet);
        bindOrder->slot = 3;
    }

    {
        auto texturesSet = bindTransmittance();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_transmittance_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        transmittanceCommands->addChild(bindPipeline);
        transmittanceCommands->addChild(bindTextures);
        transmittanceCommands->addChild(vsg::Dispatch::create(uint32_t(ceil(float(_settings->transmittanceWidth) / float(_settings->numThreads))),
                                                              uint32_t(ceil(float(_settings->transmittanceHeight) / float(_settings->numThreads))), 1));
        transmittanceCommands->addChild(memoryBarrierCmd);
        singlePass->addChild(transmittanceCommands);
    }

    {
        auto texturesSet = bindDirectIrradiance();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_direct_irradiance_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        singlePass->addChild(bindPipeline);
        singlePass->addChild(bindTextures);
        singlePass->addChild(vsg::Dispatch::create(uint32_t(ceil(float(_settings->irradianceWidth) / float(_settings->numThreads))),
                                                   uint32_t(ceil(float(_settings->irradianceHeight) / float(_settings->numThreads))), 1));
        singlePass->addChild(memoryBarrierCmd);
    }

    {
        auto texturesSet = bindSingleScattering();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_single_scattering_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        singlePass->addChild(bindPipeline);
        singlePass->addChild(bindTextures);
        singlePass->addChild(vsg::Dispatch::create(uint32_t(scatteringWidth() / _settings->numThreads),
                                                   uint32_t(scatteringHeight() / _settings->numThreads),
                                                   uint32_t(scatteringDepth() / _settings->numThreads)));
        singlePass->addChild(memoryBarrierCmd);
    }

    {
        auto texturesSet = bindScatteringDensity();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout, oLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_scattering_density_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        multipleScattering->addChild(bindPipeline);
        multipleScattering->addChild(bindTextures);
        multipleScattering->addChild(vsg::Dispatch::create(uint32_t(ceil(float(scatteringWidth()) / float(_settings->numThreads))),
                                                           uint32_t(ceil(float(scatteringHeight()) / float(_settings->numThreads))),
                                                           uint32_t(ceil(float(scatteringDepth()) / float(_settings->numThreads)))));
        multipleScattering->addChild(memoryBarrierCmd);
    }

    {
        auto texturesSet = bindIndirectIrradiance();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout, oLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_indirect_irradiance_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        multipleScattering->addChild(bindPipeline);
        multipleScattering->addChild(bindTextures);
        multipleScattering->addChild(vsg::Dispatch::create(uint32_t(ceil(float(_settings->irradianceWidth) / float(_settings->numThreads))),
                                                           uint32_t(ceil(float(_settings->irradianceHeight) / float(_settings->numThreads))), 1));
        multipleScattering->addChild(memoryBarrierCmd);
    }

    {
        auto texturesSet = bindMultipleScattering();
        auto pipelineLayout = vsg::PipelineLayout::create(vsg::DescriptorSetLayouts{texturesSet->setLayout, pLayout, oLayout}, vsg::PushConstantRanges{});
        auto bindPipeline = bindCompute("shaders/scattering/compute_multiple_scattering_cs.glsl", pipelineLayout);
        auto bindTextures = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, texturesSet);
        multipleScattering->addChild(bindPipeline);
        multipleScattering->addChild(bindTextures);
        multipleScattering->addChild(vsg::Dispatch::create(uint32_t(ceil(float(scatteringWidth()) / float(_settings->numThreads))),
                                                           uint32_t(ceil(float(scatteringHeight()) / float(_settings->numThreads))),
                                                           uint32_t(ceil(float(scatteringDepth()) / float(_settings->numThreads)))));
        multipleScattering->addChild(memoryBarrierCmd);
    }

    // compile the Vulkan objects
    auto compileTraversal = vsg::CompileTraversal::create(_device);
    auto context = compileTraversal->contexts.front();

    bindParameters->accept(*compileTraversal);
    bindOrder->accept(*compileTraversal);
    singlePass->accept(*compileTraversal);
    multipleScattering->accept(*compileTraversal);

    int computeQueueFamily = _physicalDevice->getQueueFamily(VK_QUEUE_COMPUTE_BIT);

    context->commandPool = vsg::CommandPool::create(_device, computeQueueFamily);

    auto fence = vsg::Fence::create(_device);
    auto computeQueue = _device->getQueue(computeQueueFamily);

    int num_iterations = (_settings->precomputedWavelenghts + 2) / 3;
    double dlambda = (kLambdaMax - kLambdaMin) / (3 * num_iterations);
    for (int i = 0; i < num_iterations; ++i)
    {
        vsg::vec3 lambdas(
            kLambdaMin + (3 * i + 0.5) * dlambda,
            kLambdaMin + (3 * i + 1.5) * dlambda,
            kLambdaMin + (3 * i + 2.5) * dlambda);

        auto coeff = [dlambda](double lambda, int component)
        {
            // Note that we don't include MAX_LUMINOUS_EFFICACY here, to avoid
            // artefacts due to too large values when using half precision on GPU.
            // We add this term back in kAtmosphereShader, via
            // SKY_SPECTRAL_RADIANCE_TO_LUMINANCE (see also the comments in the
            // Model constructor).
            double x = cieColorMatchingFunctionTableValue(lambda, 1);
            double y = cieColorMatchingFunctionTableValue(lambda, 2);
            double z = cieColorMatchingFunctionTableValue(lambda, 3);
            return static_cast<float>((
                                          XYZ_TO_SRGB[component * 3] * x +
                                          XYZ_TO_SRGB[component * 3 + 1] * y +
                                          XYZ_TO_SRGB[component * 3 + 2] * z) * dlambda);
        };
        vsg::mat4 luminance_from_radiance{
                                          coeff(lambdas[0], 0), coeff(lambdas[0], 1), coeff(lambdas[0], 2), 0.0f,
                                          coeff(lambdas[1], 0), coeff(lambdas[1], 1), coeff(lambdas[1], 2), 0.0f,
                                          coeff(lambdas[2], 0), coeff(lambdas[2], 1), coeff(lambdas[2], 2), 0.0f,
                                          0.0, 0.0f, 0.0f, 1.0f};

        lfr->set(luminance_from_radiance);
        computeParameters(parameters, lambdas);

        lfrBuffer->copyDataListToBuffers();
        parametersBuffer->copyDataListToBuffers();

        vsg::submitCommandsToQueue(context->commandPool, fence, 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
            bindParameters->record(commandBuffer);
            singlePass->record(commandBuffer);
        });

        // Compute the 2nd, 3rd and 4th orderValue of scattering, in sequence.
        for (int o = 2; o <= _settings->scatteringOrders; ++o)
        {
            order->set(o);
            orderBuffer->copyDataListToBuffers();

            vsg::submitCommandsToQueue(context->commandPool, fence, 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
                bindParameters->record(commandBuffer);
                bindOrder->record(commandBuffer);
                multipleScattering->record(commandBuffer);
            });
        }
    }

    lfr->set(vsg::mat4());
    computeParameters(parameters, {kLambdaR, kLambdaG, kLambdaB});

    lfrBuffer->copyDataListToBuffers();
    parametersBuffer->copyDataListToBuffers();

    vsg::submitCommandsToQueue(context->commandPool, fence, 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
        bindParameters->record(commandBuffer);
        transmittanceCommands->record(commandBuffer);
    });
}

vsg::vec3 AtmosphereGenerator::convertSpectrumToLinearSrgb(double c)
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    const int dlambda = 1;
    for (int lambda = kLambdaMin; lambda < kLambdaMax; lambda += dlambda) {
      double value = interpolate(waveLengths, solarIrradiance, lambda);
      x += cieColorMatchingFunctionTableValue(lambda, 1) * value;
      y += cieColorMatchingFunctionTableValue(lambda, 2) * value;
      z += cieColorMatchingFunctionTableValue(lambda, 3) * value;
    }
    auto r = MAX_LUMINOUS_EFFICACY *
        (XYZ_TO_SRGB[0] * x + XYZ_TO_SRGB[1] * y + XYZ_TO_SRGB[2] * z) * dlambda;
    auto g = MAX_LUMINOUS_EFFICACY *
        (XYZ_TO_SRGB[3] * x + XYZ_TO_SRGB[4] * y + XYZ_TO_SRGB[5] * z) * dlambda;
    auto b = MAX_LUMINOUS_EFFICACY *
        (XYZ_TO_SRGB[6] * x + XYZ_TO_SRGB[7] * y + XYZ_TO_SRGB[8] * z) * dlambda;

    double white_point = (r + g + b) / c;

    r /= white_point;
    g /= white_point;
    b /= white_point;

    return {static_cast<float>(r), static_cast<float>(g), static_cast<float>(b)};
}

vsg::ref_ptr<AtmosphereBinding> AtmosphereGenerator::loadData()
{
    auto atmosphereBinding = AtmosphereBinding::create();
    atmosphereBinding->transmittanceTexture = _transmittanceTexture;
    atmosphereBinding->irradianceTexture = _irradianceTexture;
    atmosphereBinding->scatteringTexture = _scatteringTexture;
    atmosphereBinding->singleMieScatteringTexture = _singleMieScatteringTexture;

    return atmosphereBinding;
}

vsg::ref_ptr<AtmosphereRuntime> AtmosphereGenerator::createRuntime(vsg::ref_ptr<AtmosphereBinding> atmosphere, vsg::ref_ptr<CloudsBinding> clouds)
{
    auto runtimeData = AtmosphereRuntime::create(atmosphere, clouds);
    runtimeData->createPhongShaderSet(_options, _renderConstants);
    runtimeData->createSkyShaderSet(_options, _renderConstants);

    runtimeData->ellipsoidModel = _settings->ellipsoidModel;
    runtimeData->atmosphereBinding->settings->value() = {vsg::vec4(convertSpectrumToLinearSrgb(3.0), 0.0f), {std::tan(_settings->sunAngularRadius), std::cos(_settings->sunAngularRadius)}};
    runtimeData->lengthUnitInMeters = _settings->lengthUnitInMeters;

    return runtimeData;
}

void AtmosphereGenerator::generateTextures()
{
    auto sampler = vsg::Sampler::create();
    _transmittanceTexture = Image::create(VkExtent3D{_settings->transmittanceWidth, _settings->transmittanceHeight, 1}, sampler);
    _transmittanceTexture->allocateTexture(_device);

    _irradianceTexture = Image::create(VkExtent3D{_settings->irradianceWidth, _settings->irradianceHeight, 1}, sampler);
    _irradianceTexture->allocateTexture(_device, true);
    _deltaIrradianceTexture = Image::create(VkExtent3D{_settings->irradianceWidth, _settings->irradianceHeight, 1}, sampler);
    _deltaIrradianceTexture->allocateTexture(_device);

    VkExtent3D extent{scatteringWidth(), scatteringHeight(), scatteringDepth()};

    _deltaRayleighScatteringTexture = Image::create(extent, sampler);
    _deltaRayleighScatteringTexture->allocateTexture(_device);
    _deltaMieScatteringTexture = Image::create(extent, sampler);
    _deltaMieScatteringTexture->allocateTexture(_device);
    _scatteringTexture = Image::create(extent, sampler);
    _scatteringTexture->allocateTexture(_device, true);
    _singleMieScatteringTexture = Image::create(extent, sampler);
    _singleMieScatteringTexture->allocateTexture(_device, true);

    _deltaScatteringDensityTexture = Image::create(extent, sampler);
    _deltaScatteringDensityTexture->allocateTexture(_device);
    _deltaMultipleScatteringTexture = Image::create(extent, sampler);
    _deltaMultipleScatteringTexture->allocateTexture(_device);
}

vsg::ref_ptr<vsg::BindComputePipeline> AtmosphereGenerator::bindCompute(const vsg::Path& filename, vsg::ref_ptr<vsg::PipelineLayout> pipelineLayout) const
{
    auto computeStage = vsg::ShaderStage::read(VK_SHADER_STAGE_COMPUTE_BIT, "main", filename, _options);
    if (!computeStage)
    {
        vsg::error("Could not find shaders.");
        return {};
    }

    computeStage->module->hints = compileSettings;
    computeStage->specializationConstants = _computeConstants;

    // set up the compute pipeline
    auto pipeline = vsg::ComputePipeline::create(pipelineLayout, computeStage);

    return vsg::BindComputePipeline::create(pipeline);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindTransmittance() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{{0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}};
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    auto writeTexture = vsg::DescriptorImage::create(_transmittanceTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    vsg::Descriptors descriptors{writeTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindDirectIrradiance() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
                                                        {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}};
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    auto deltaTexture = vsg::DescriptorImage::create(_deltaIrradianceTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    auto transmittanceTexture = vsg::DescriptorImage::create(_transmittanceTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);

    vsg::Descriptors descriptors{deltaTexture, transmittanceTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindSingleScattering() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
                                                        {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {1, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {2, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {3, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}};
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    auto rayTexture = vsg::DescriptorImage::create(_deltaRayleighScatteringTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    auto mieTexture = vsg::DescriptorImage::create(_deltaMieScatteringTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    auto readTexture = vsg::DescriptorImage::create(_scatteringTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    auto writeTexture = vsg::DescriptorImage::create(_singleMieScatteringTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    auto transmittanceTexture = vsg::DescriptorImage::create(_transmittanceTexture->imageInfo, 4, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);

    vsg::Descriptors descriptors{rayTexture, mieTexture, readTexture, writeTexture, transmittanceTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindScatteringDensity() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
        {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {5, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}
    };
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    auto scatteringTexture = vsg::DescriptorImage::create(_deltaScatteringDensityTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    auto transmittanceTexture = vsg::DescriptorImage::create(_transmittanceTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto rayTexture = vsg::DescriptorImage::create(_deltaRayleighScatteringTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto mieTexture = vsg::DescriptorImage::create(_deltaMieScatteringTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto multipleTexture = vsg::DescriptorImage::create(_deltaMultipleScatteringTexture->imageInfo, 4, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto irradianceTexture = vsg::DescriptorImage::create(_deltaIrradianceTexture->imageInfo, 5, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);

    vsg::Descriptors descriptors{scatteringTexture, transmittanceTexture, rayTexture, mieTexture, multipleTexture, irradianceTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindIndirectIrradiance() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
                                                        {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {1, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        };
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

    auto deltaIrradianceTexture = vsg::DescriptorImage::create(_deltaIrradianceTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    auto irradianceTexture = vsg::DescriptorImage::create(_irradianceTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    auto singleTexture = vsg::DescriptorImage::create(_deltaRayleighScatteringTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto mieTexture = vsg::DescriptorImage::create(_deltaMieScatteringTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto multipleTexture = vsg::DescriptorImage::create(_deltaMultipleScatteringTexture->imageInfo, 4, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);

    //auto orderValueBuffer = vsg::DescriptorBuffer::create(orderValue, 6);

    vsg::Descriptors descriptors{deltaIrradianceTexture, irradianceTexture, singleTexture,
                                 mieTexture, multipleTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSet> AtmosphereGenerator::bindMultipleScattering() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
        {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {1, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
        {3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}
    };
    auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);


    auto deltaMultipleTexture = vsg::DescriptorImage::create(_deltaMultipleScatteringTexture->imageInfo, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);
    auto scatteringTexture = vsg::DescriptorImage::create(_scatteringTexture->imageInfo, 1, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

    auto transmittanceTexture = vsg::DescriptorImage::create(_transmittanceTexture->imageInfo, 2, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);
    auto deltaDensityTexture = vsg::DescriptorImage::create(_deltaScatteringDensityTexture->imageInfo, 3, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER);

    vsg::Descriptors descriptors{deltaMultipleTexture, scatteringTexture, transmittanceTexture, deltaDensityTexture};
    return vsg::DescriptorSet::create(descriptorSetLayout, descriptors);
}

vsg::ref_ptr<vsg::DescriptorSetLayout> AtmosphereGenerator::parametersLayout() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{
                                                        {0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        {2, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr},
                                                        };
    return vsg::DescriptorSetLayout::create(descriptorBindings);
}

vsg::ref_ptr<vsg::DescriptorSetLayout> AtmosphereGenerator::orderLayout() const
{
    // set up DescriptorSetLayout, DecriptorSet and BindDescriptorSets
    vsg::DescriptorSetLayoutBindings descriptorBindings{{0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}};
    return vsg::DescriptorSetLayout::create(descriptorBindings);
}

void AtmosphereGenerator::computeParameters(vsg::ref_ptr<vsg::Value<Parameters>> parameters, const vsg::vec3 &lambdas) const
{
    Parameters p;
    p.solar_irradiance = toVector(waveLengths, solarIrradiance, lambdas, 1.0);
    p.rayleigh_scattering = toVector(waveLengths, rayleighScattering, lambdas, _settings->lengthUnitInMeters);
    p.mie_scattering = toVector(waveLengths, mieScattering, lambdas, _settings->lengthUnitInMeters);
    p.mie_extinction = toVector(waveLengths, mieExtinction, lambdas, _settings->lengthUnitInMeters);
    p.absorption_extinction = toVector(waveLengths, absorptionExtinction, lambdas, _settings->lengthUnitInMeters);
    p.ground_albedo = toVector(waveLengths, groundAlbedo, lambdas, 1.0);

    parameters->set(p);
}

// -----------------------------------------------------------------------------------------------------------------------------------

void AtmosphereGenerator::assignRenderConstants()
{
    auto SKY_SPECTRAL_RADIANCE_TO_LUMINANCE = computeSpectralRadianceToLuminanceFactors(waveLengths, solarIrradiance, -3);
    auto SUN_SPECTRAL_RADIANCE_TO_LUMINANCE = computeSpectralRadianceToLuminanceFactors(waveLengths, solarIrradiance, 0);

    vsg::vec3 lambdas(kLambdaR, kLambdaG, kLambdaB);

    auto solar_irradiance = toVector(waveLengths, solarIrradiance, lambdas, 1.0);
    auto rayleigh_scattering = toVector(waveLengths, rayleighScattering, lambdas, _settings->lengthUnitInMeters);
    auto mie_scattering = toVector(waveLengths, mieScattering, lambdas, _settings->lengthUnitInMeters);
    auto mie_extinction = toVector(waveLengths, mieExtinction, lambdas, _settings->lengthUnitInMeters);
    auto absorption_extinction = toVector(waveLengths, absorptionExtinction, lambdas, _settings->lengthUnitInMeters);
    auto ground_albedo = toVector(waveLengths, groundAlbedo, lambdas, 1.0);

    _renderConstants = {

        {0, vsg::intValue::create(_settings->numThreads)},
        //{1, vsg::intValue::create(_settings->numViewerThreads)},
        //{2, vsg::intValue::create(_settings->cubeSize)},

        {10, vsg::intValue::create(_settings->transmittanceWidth)},
        {11, vsg::intValue::create(_settings->transmittanceHeight)},

        {12, vsg::intValue::create(_settings->scaterringR)},
        {13, vsg::intValue::create(_settings->scaterringMU)},
        {14, vsg::intValue::create(_settings->scaterringMU_S)},
        {15, vsg::intValue::create(_settings->scaterringNU)},

        {16, vsg::intValue::create(_settings->irradianceWidth)},
        {17, vsg::intValue::create(_settings->irradianceHeight)},

        {18, vsg::floatValue::create(rayleigh_scattering.x)},
        {19, vsg::floatValue::create(rayleigh_scattering.y)},
        {20, vsg::floatValue::create(rayleigh_scattering.z)},

        {21, vsg::floatValue::create(mie_scattering.x)},
        {22, vsg::floatValue::create(mie_scattering.y)},
        {23, vsg::floatValue::create(mie_scattering.z)},

        {24, vsg::floatValue::create(absorption_extinction.x)},
        {25, vsg::floatValue::create(absorption_extinction.y)},
        {26, vsg::floatValue::create(absorption_extinction.z)},

        {27, vsg::floatValue::create(solar_irradiance.x)},
        {28, vsg::floatValue::create(solar_irradiance.y)},
        {29, vsg::floatValue::create(solar_irradiance.z)},

        {30, vsg::floatValue::create(mie_extinction.x)},
        {31, vsg::floatValue::create(mie_extinction.y)},
        {32, vsg::floatValue::create(mie_extinction.z)},

        {33, vsg::floatValue::create(ground_albedo.x)},
        {34, vsg::floatValue::create(ground_albedo.y)},
        {35, vsg::floatValue::create(ground_albedo.z)},

        {36, vsg::floatValue::create(_settings->sunAngularRadius)},
        {37, vsg::floatValue::create(static_cast<float>(_settings->ellipsoidModel->radiusEquator() / _settings->lengthUnitInMeters))},
        {38, vsg::floatValue::create(static_cast<float>((_settings->ellipsoidModel->radiusEquator() + _settings->atmoshpereHeight) / _settings->lengthUnitInMeters))},
        {39, vsg::floatValue::create(_settings->miePhaseFunction_g)},
        {40, vsg::floatValue::create(static_cast<float>(std::cos(_settings->maxSunZenithAngle)))},

        {41, vsg::floatValue::create(SKY_SPECTRAL_RADIANCE_TO_LUMINANCE.x)},
        {42, vsg::floatValue::create(SKY_SPECTRAL_RADIANCE_TO_LUMINANCE.y)},
        {43, vsg::floatValue::create(SKY_SPECTRAL_RADIANCE_TO_LUMINANCE.z)},

        {44, vsg::floatValue::create(SUN_SPECTRAL_RADIANCE_TO_LUMINANCE.x)},
        {45, vsg::floatValue::create(SUN_SPECTRAL_RADIANCE_TO_LUMINANCE.y)},
        {46, vsg::floatValue::create(SUN_SPECTRAL_RADIANCE_TO_LUMINANCE.z)}
    };

}

vsg::ref_ptr<AtmosphereGenerator> createAtmosphereGenerator(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::EllipsoidModel> eps, vsg::ref_ptr<vsg::Options> options)
{
    // Values from "Reference Solar Spectral Irradiance: ASTM G-173", ETR column
    // (see http://rredc.nrel.gov/solar/spectra/am1.5/ASTMG173/ASTMG173.html),
    // summed and averaged in each bin (e.g. the value for 360nm is the average
    // of the ASTM G-173 values for all wavelengths between 360 and 370nm).
    // Values in W.m^-2.
    constexpr int kLambdaMinN = static_cast<int>(kLambdaMin);
    constexpr int kLambdaMaxN = static_cast<int>(kLambdaMax);;
    constexpr double kSolarIrradiance[48] = {
        1.11776, 1.14259, 1.01249, 1.14716, 1.72765, 1.73054, 1.6887, 1.61253,
        1.91198, 2.03474, 2.02042, 2.02212, 1.93377, 1.95809, 1.91686, 1.8298,
        1.8685, 1.8931, 1.85149, 1.8504, 1.8341, 1.8345, 1.8147, 1.78158, 1.7533,
        1.6965, 1.68194, 1.64654, 1.6048, 1.52143, 1.55622, 1.5113, 1.474, 1.4482,
        1.41018, 1.36775, 1.34188, 1.31429, 1.28303, 1.26758, 1.2367, 1.2082,
        1.18737, 1.14683, 1.12362, 1.1058, 1.07124, 1.04992
    };
    // Values from http://www.iup.uni-bremen.de/gruppen/molspec/databases/
    // referencespectra/o3spectra2011/index.html for 233K, summed and averaged in
    // each bin (e.g. the value for 360nm is the average of the original values
    // for all wavelengths between 360 and 370nm). Values in m^2.
    constexpr double kOzoneCrossSection[48] = {
        1.18e-27, 2.182e-28, 2.818e-28, 6.636e-28, 1.527e-27, 2.763e-27, 5.52e-27,
        8.451e-27, 1.582e-26, 2.316e-26, 3.669e-26, 4.924e-26, 7.752e-26, 9.016e-26,
        1.48e-25, 1.602e-25, 2.139e-25, 2.755e-25, 3.091e-25, 3.5e-25, 4.266e-25,
        4.672e-25, 4.398e-25, 4.701e-25, 5.019e-25, 4.305e-25, 3.74e-25, 3.215e-25,
        2.662e-25, 2.238e-25, 1.852e-25, 1.473e-25, 1.209e-25, 9.423e-26, 7.455e-26,
        6.566e-26, 5.105e-26, 4.15e-26, 4.228e-26, 3.237e-26, 2.451e-26, 2.801e-26,
        2.534e-26, 1.624e-26, 1.465e-26, 2.078e-26, 1.383e-26, 7.105e-27
    };
    // From https://en.wikipedia.org/wiki/Dobson_unit, in molecules.m^-2.
    constexpr double kDobsonUnit = 2.687e20;
    // Maximum number density of ozone molecules, in m^-3 (computed so at to get
    // 300 Dobson units of ozone - for this we divide 300 DU by the integral of
    // the ozone density profile defined below, which is equal to 15km).
    constexpr double kMaxOzoneNumberDensity = 300.0 * kDobsonUnit / 15000.0;
    // Wavelength independent solar irradiance "spectrum" (not physically
    // realistic, but was used in the original implementation).
    constexpr double kConstantSolarIrradiance = 1.5;
    constexpr double kRayleigh = 1.24062e-6;
    constexpr double kRayleighScaleHeight = 8000.0;
    constexpr double kMieScaleHeight = 1200.0;
    constexpr double kMieAngstromAlpha = 0.0;
    constexpr double kMieAngstromBeta = 5.328e-3;
    constexpr double kMieSingleScatteringAlbedo = 0.9;
    constexpr float kMiePhaseFunctionG = 0.8f;
    constexpr double kGroundAlbedo = 0.1;
    const double max_sun_zenith_angle = 120.0 / 180.0 * vsg::PI;

    atmosphere::DensityProfileLayer rayleigh_layer(0.0, 1.0, -1.0 / kRayleighScaleHeight, 0.0, 0.0, 1000.0);
    atmosphere::DensityProfileLayer mie_layer(0.0f, 1.0f, -1.0f / kMieScaleHeight, 0.0f, 0.0f, 1000.0);

    atmosphere::DensityProfileLayer absorption_layer0(25000.0, 0.0, 0.0, 1.0 / 15000.0, -2.0 / 3.0, 1000.0);
    atmosphere::DensityProfileLayer absorption_layer1(0.0, 0.0, 0.0, -1.0 / 15000.0, 8.0 / 3.0, 1000.0);

    // Density profile increasing linearly from 0 to 1 between 10 and 25km, and
    // decreasing linearly from 1 to 0 between 25 and 40km. This is an approximate
    // profile from http://www.kln.ac.lk/science/Chemistry/Teaching_Resources/
    // Documents/Introduction%20to%20atmospheric%20chemistry.pdf (page 10).

    bool use_constant_solar_spectrum_ = false;
    bool use_ozone = true;

    std::vector<double> wavelengths;
    std::vector<double> solar_irradiance;
    std::vector<double> rayleigh_scattering;
    std::vector<double> mie_scattering;
    std::vector<double> mie_extinction;
    std::vector<double> absorption_extinction;
    std::vector<double> ground_albedo;
    for (int l = kLambdaMinN; l <= kLambdaMaxN; l += 10) {
        double lambda = static_cast<double>(l) * 1e-3;  // micro-meters
        double mie =
                kMieAngstromBeta / kMieScaleHeight * pow(lambda, -kMieAngstromAlpha);
        wavelengths.push_back(l);
        if (use_constant_solar_spectrum_) {
            solar_irradiance.push_back(kConstantSolarIrradiance);
        } else {
            solar_irradiance.push_back(kSolarIrradiance[(l - kLambdaMinN) / 10]);
        }
        rayleigh_scattering.push_back(kRayleigh * pow(lambda, -4));
        mie_scattering.push_back(mie * kMieSingleScatteringAlbedo);
        mie_extinction.push_back(mie);
        absorption_extinction.push_back(use_ozone ?
                                            kMaxOzoneNumberDensity * kOzoneCrossSection[(l - kLambdaMinN) / 10] :
                                        0.0);
        ground_albedo.push_back(kGroundAlbedo);
    }

    auto model = atmosphere::AtmosphereGenerator::create(window->getOrCreateDevice(), window->getOrCreatePhysicalDevice(), options);
/*
    model->waveLengths = wavelengths;
    model->solarIrradiance = solar_irradiance;
    model->sunAngularRadius = 0.01935f;
    model->ellipsoidModel = eps;
    model->atmoshpereHeight = 60000.0;        auto atmosphere = atmosphereGenerator->loadData();
    model->rayleighDensityLayer = rayleigh_layer;
    model->rayleighScattering = rayleigh_scattering;
    model->mieDensityLayer = mie_layer;
    model->mieScattering = mie_scattering;
    model->mieExtinction = mie_extinction;
    model->miePhaseFunction_g = kMiePhaseFunctionG;
    model->absorptionDensityLayer0 = absorption_layer0;
    model->absorptionDensityLayer1 = absorption_layer1;
    model->absorptionExtinction = absorption_extinction;
    model->groundAlbedo = ground_albedo;
    model->maxSunZenithAngle = max_sun_zenith_angle;
    model->lengthUnitInMeters = 1000.0;

    model->initialize(4);
*/
    return model;
}

vsg::ref_ptr<AtmosphereGenerator> createAtmosphereGenerator(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Options> options)
{
    auto model = atmosphere::AtmosphereGenerator::create(settings, window->getOrCreateDevice(), window->getOrCreatePhysicalDevice(), options);

    model->initialize();

    return model;
}

vsg::ref_ptr<AtmosphereGenerator> createAtmosphereGenerator(vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Options> options)
{
    vsg::Names instanceExtensions;
    vsg::Names requestedLayers;
    vsg::Names deviceExtensions;

    vsg::Names validatedNames = vsg::validateInstancelayerNames(requestedLayers);

    // get the physical device that supports the required compute queue
    auto instance = vsg::Instance::create(instanceExtensions, validatedNames);
    auto [physicalDevice, computeQueueFamily] = instance->getPhysicalDeviceAndQueueFamily(VK_QUEUE_COMPUTE_BIT);
    if (!physicalDevice || computeQueueFamily < 0)
    {
        vsg::error("No vkPhysicalDevice available that supports compute.");
        return {};
    }

    // create the logical device with specified queue, layers and extensions
    vsg::QueueSettings queueSettings{vsg::QueueSetting{computeQueueFamily, {1.0}}};
    auto device = vsg::Device::create(physicalDevice, queueSettings, validatedNames, deviceExtensions);

    auto model = atmosphere::AtmosphereGenerator::create(settings, device, physicalDevice, options);

    model->initialize();

    return model;
}

vsg::ref_ptr<AtmosphereRuntime> createAtmosphereRuntime(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<AtmosphereModelSettings> settings, vsg::ref_ptr<vsg::Options> options)
{
    auto model = atmosphere::AtmosphereGenerator::create(settings, window->getOrCreateDevice(), window->getOrCreatePhysicalDevice(), options);
    model->initialize();
    auto atmosphere = model->loadData();

    vsg::ref_ptr<CloudsBinding> clouds;
    if(settings->clouds)
        clouds = createCloudsData(window, options, settings);

    return model->createRuntime(atmosphere, clouds);
}

}
