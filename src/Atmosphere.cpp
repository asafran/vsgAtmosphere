#include "Atmosphere.h"
#include "AtmosphereTools.h"
#include "BRDF.h"
#include "Clouds.h"
#include "Task.h"

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

    _computeConstants = _settings->getComputeConstants();
    assignRenderConstants();

    auto memoryBarrier = vsg::MemoryBarrier::create();
    memoryBarrier->srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    memoryBarrier->dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
    auto memoryBarrierCmd = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, memoryBarrier);

    auto parameters = ComputeParametersBinding::create(_settings);
    auto singlePass = vsg::Commands::create();
    auto multipleScattering = vsg::Commands::create();

    auto transitionBarrierCmd = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0);
    transitionBarrierCmd->add(_transmittanceTexture->transitionWriteBarrier());

    transitionBarrierCmd->add(_irradianceTexture->transitionWriteBarrier());
    transitionBarrierCmd->add(_deltaIrradianceTexture->transitionWriteBarrier());

    transitionBarrierCmd->add(_deltaRayleighScatteringTexture->transitionWriteBarrier());
    transitionBarrierCmd->add(_deltaMieScatteringTexture->transitionWriteBarrier());
    transitionBarrierCmd->add(_scatteringTexture->transitionWriteBarrier());
    transitionBarrierCmd->add(_singleMieScatteringTexture->transitionWriteBarrier());

    transitionBarrierCmd->add(_deltaScatteringDensityTexture->transitionWriteBarrier());
    transitionBarrierCmd->add(_deltaMultipleScatteringTexture->transitionWriteBarrier());

    auto transmittanceTask = Task::create();
    transmittanceTask->shader = createComputeShader(std::string(transmittanceShader));
    transmittanceTask->numThreads = _settings->numThreads;
    transmittanceTask->parameters = parameters;
    transmittanceTask->writeImages = {_transmittanceTexture};
    transmittanceTask->readImages = {};

    auto directIrradianceTask = Task::create();
    directIrradianceTask->shader = createComputeShader(std::string(directIrradianceShader));
    directIrradianceTask->numThreads = _settings->numThreads;
    directIrradianceTask->parameters = parameters;
    directIrradianceTask->writeImages = {_deltaIrradianceTexture};
    directIrradianceTask->readImages = {_transmittanceTexture};

    auto singleScatteringTask = Task::create();
    singleScatteringTask->shader = createComputeShader(std::string(singleScatteringShader));
    singleScatteringTask->numThreads = _settings->numThreads;
    singleScatteringTask->parameters = parameters;
    singleScatteringTask->writeImages = {_deltaRayleighScatteringTexture, _deltaMieScatteringTexture, _scatteringTexture, _singleMieScatteringTexture};
    singleScatteringTask->readImages = {_transmittanceTexture};

    auto scatteringDensityTask = Task::create();
    scatteringDensityTask->shader = createComputeShader(std::string(scatteringDensityShader));
    scatteringDensityTask->numThreads = _settings->numThreads;
    scatteringDensityTask->parameters = parameters;
    scatteringDensityTask->writeImages = {_deltaScatteringDensityTexture};
    scatteringDensityTask->readImages = {_transmittanceTexture,
                                         _deltaRayleighScatteringTexture,
                                         _deltaMieScatteringTexture,
                                         _deltaMultipleScatteringTexture,
                                         _deltaIrradianceTexture};

    auto indirectIrradianceTask = Task::create();
    indirectIrradianceTask->shader = createComputeShader(std::string(indirectIrradianceShader));
    indirectIrradianceTask->numThreads = _settings->numThreads;
    indirectIrradianceTask->parameters = parameters;
    indirectIrradianceTask->writeImages = {_deltaIrradianceTexture, _irradianceTexture};
    indirectIrradianceTask->readImages = {_deltaRayleighScatteringTexture, _deltaMieScatteringTexture, _deltaMultipleScatteringTexture};

    auto multipleScatteringTask = Task::create();
    multipleScatteringTask->shader = createComputeShader(std::string(multipleScatteringShader));
    multipleScatteringTask->numThreads = _settings->numThreads;
    multipleScatteringTask->parameters = parameters;
    multipleScatteringTask->writeImages = {_deltaMultipleScatteringTexture, _scatteringTexture};
    multipleScatteringTask->readImages = {_transmittanceTexture, _deltaScatteringDensityTexture};

    singlePass->addChild(transmittanceTask->createTaskCommands());
    singlePass->addChild(memoryBarrierCmd);
    singlePass->addChild(directIrradianceTask->createTaskCommands());
    singlePass->addChild(memoryBarrierCmd);
    singlePass->addChild(singleScatteringTask->createTaskCommands());
    singlePass->addChild(memoryBarrierCmd);

    multipleScattering->addChild(scatteringDensityTask->createTaskCommands());
    multipleScattering->addChild(memoryBarrierCmd);
    multipleScattering->addChild(indirectIrradianceTask->createTaskCommands());
    multipleScattering->addChild(memoryBarrierCmd);
    multipleScattering->addChild(multipleScatteringTask->createTaskCommands());
    multipleScattering->addChild(memoryBarrierCmd);

    for (int o = 2; o <= _settings->scatteringOrders; ++o)
    {
        auto order = vsg::intValue::create(o);
        auto pushOrder = vsg::PushConstants::create(VK_SHADER_STAGE_COMPUTE_BIT, 0, order);
        singlePass->addChild(pushOrder);
        singlePass->addChild(multipleScattering);
    }

    // compile the Vulkan objects
    auto compileTraversal = vsg::CompileTraversal::create(_device);
    auto context = compileTraversal->contexts.front();

    transitionBarrierCmd->accept(*compileTraversal);
    singlePass->accept(*compileTraversal);

    int computeQueueFamily = _physicalDevice->getQueueFamily(VK_QUEUE_COMPUTE_BIT);
    context->commandPool = vsg::CommandPool::create(_device, computeQueueFamily);
    auto computeQueue = _device->getQueue(computeQueueFamily);

    vsg::submitCommandsToQueue(context->commandPool, vsg::Fence::create(_device), 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
        transitionBarrierCmd->record(commandBuffer);
    });

    if(_settings->radiance)
    {
        int num_iterations = (_settings->precomputedWavelenghts + 2) / 3;
        double dlambda = (kLambdaMax - kLambdaMin) / (3 * num_iterations);
        for (int i = 0; i < num_iterations; ++i)
        {
            vsg::vec3 lambdas(
                kLambdaMin + (3 * i + 0.5) * dlambda,
                kLambdaMin + (3 * i + 1.5) * dlambda,
                kLambdaMin + (3 * i + 2.5) * dlambda);

            parameters->computeLuminanceFromRadiance(lambdas, dlambda);
            parameters->setParameters(computeParameters(lambdas));

            vsg::submitCommandsToQueue(context->commandPool, vsg::Fence::create(_device), 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
                singlePass->record(commandBuffer);
            });
        }
    }

    parameters->resetLuminanceFromRadiance();
    parameters->setParameters(computeParameters({kLambdaR, kLambdaG, kLambdaB}));

    if(!_settings->radiance)
    {
        vsg::submitCommandsToQueue(context->commandPool, vsg::Fence::create(_device), 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
            singlePass->record(commandBuffer);
            multipleScattering->record(commandBuffer);
        });
    }
    else
    {
        vsg::submitCommandsToQueue(context->commandPool, vsg::Fence::create(_device), 100000000000, computeQueue, [&](vsg::CommandBuffer& commandBuffer) {
            transmittanceTask->createTaskCommands()->record(commandBuffer);
        });
    }

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

    atmosphereBinding->radiance = _settings->radiance;

    return atmosphereBinding;
}

vsg::ref_ptr<AtmosphereRuntime> AtmosphereGenerator::createRuntime(vsg::ref_ptr<AtmosphereBinding> atmosphere, vsg::ref_ptr<BRDFBinding> pbr, vsg::ref_ptr<CloudsBinding> clouds)
{
    auto runtimeData = AtmosphereRuntime::create(atmosphere, pbr, clouds);
    runtimeData->createPhongShaderSet(_options, _renderConstants);
    runtimeData->createSkyShaderSet(_options, _renderConstants);

    runtimeData->ellipsoidModel = _settings->ellipsoidModel;
    runtimeData->atmosphereBinding->settings->value() = {vsg::vec4(convertSpectrumToLinearSrgb(3.0), 0.0f), {std::tan(_settings->sunAngularRadius), std::cos(_settings->sunAngularRadius)}};
    runtimeData->lengthUnitInMeters = _settings->lengthUnitInMeters;

    runtimeData->exposureModifier = _settings->radiance ? 1e-6 : 1.0;

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

vsg::ref_ptr<vsg::ShaderStage> AtmosphereGenerator::createComputeShader(const std::string &key) const
{
    auto shaderModule = _options->getRefObject<vsg::ShaderModule>(key);
    auto computeStage = vsg::ShaderStage::create(VK_SHADER_STAGE_COMPUTE_BIT, "main", shaderModule);
    if (!computeStage)
    {
        vsg::error("Could not find shader:", key);
        return {};
    }

    computeStage->module->hints = compileSettings;
    computeStage->specializationConstants = _computeConstants;
    return computeStage;
}

Parameters AtmosphereGenerator::computeParameters(const vsg::vec3 &lambdas) const
{
    Parameters p;
    p.solar_irradiance = toVector(waveLengths, solarIrradiance, lambdas, 1.0);
    p.rayleigh_scattering = toVector(waveLengths, rayleighScattering, lambdas, _settings->lengthUnitInMeters);
    p.mie_scattering = toVector(waveLengths, mieScattering, lambdas, _settings->lengthUnitInMeters);
    p.mie_extinction = toVector(waveLengths, mieExtinction, lambdas, _settings->lengthUnitInMeters);
    p.absorption_extinction = toVector(waveLengths, absorptionExtinction, lambdas, _settings->lengthUnitInMeters);
    p.ground_albedo = toVector(waveLengths, groundAlbedo, lambdas, 1.0);

    return p;
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
    auto atmosphereBinding = model->loadData();

    vsg::ref_ptr<CloudsBinding> cloudsBinding;
    if(settings->clouds)
        cloudsBinding = createData<CloudsGenerator, CloudsBinding>(window, options, settings);

    auto pbrBinding = createData<BRDFGenerator, BRDFBinding>(window, options, settings);

    return model->createRuntime(atmosphereBinding, pbrBinding, cloudsBinding);
}

}
