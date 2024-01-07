#include "AtmosphereRuntime.h"
#include "AtmosphereTools.h"
#include "InverseMatrices.h"

#include <vsg/all.h>

namespace atmosphere {

    AtmosphereRuntime::AtmosphereRuntime(vsg::ref_ptr<AtmosphereBinding> atmosphere, vsg::ref_ptr<BRDFBinding> pbr, vsg::ref_ptr<CloudsBinding> clouds)
        : atmosphereBinding(atmosphere)
        , cloudsBinding(clouds)
        , pbrBinding(pbr)
    {
        positionalBinding = PositionalBinding::create();
        inversePositionalBinding = PositionalBinding::create();
    }

    AtmosphereRuntime::AtmosphereRuntime()
    {

    }

    int AtmosphereRuntime::compare(const Object &rhs) const
    {
        return vsg::Object::compare(rhs);
    }

    bool AtmosphereRuntime::createPhongShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts)
    {
        auto vertexModule = options->getRefObject<vsg::ShaderModule>(std::string(phongVertShader));
        auto fragmentModule = options->getRefObject<vsg::ShaderModule>(std::string(phongFragShader));

        if (!vertexModule || !fragmentModule)
        {
            vsg::warn("Could not find shaders. Phong support disabled");
            return false;
        }
        else
            vsg::info("Created PhongShaderSet");

        auto vertexShader = vsg::ShaderStage::create(VK_SHADER_STAGE_VERTEX_BIT, "main", vertexModule);
        auto fragmentShader = vsg::ShaderStage::create(VK_SHADER_STAGE_FRAGMENT_BIT, "main", fragmentModule);
        fragmentShader->specializationConstants = constatnts;

        auto shaderSet = vsg::ShaderSet::create(vsg::ShaderStages{vertexShader, fragmentShader});

        shaderSet->addAttributeBinding("vsg_Vertex", "", 0, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_Normal", "", 1, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_TexCoord0", "", 2, VK_FORMAT_R32G32_SFLOAT, vsg::vec2Array::create(1));
        shaderSet->addAttributeBinding("vsg_Color", "", 3, VK_FORMAT_R32G32B32A32_SFLOAT, vsg::vec4Array::create(1));

        shaderSet->addAttributeBinding("vsg_position", "VSG_INSTANCE_POSITIONS", 4, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_position_scaleDistance", "VSG_BILLBOARD", 4, VK_FORMAT_R32G32B32A32_SFLOAT, vsg::vec4Array::create(1));

        shaderSet->addDescriptorBinding("displacementMap", "VSG_DISPLACEMENT_MAP", MATERIAL_DESCRIPTOR_SET, 6, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT, vsg::floatArray2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));
        shaderSet->addDescriptorBinding("diffuseMap", "VSG_DIFFUSE_MAP", MATERIAL_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::ubvec4Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R8G8B8A8_UNORM}));
        shaderSet->addDescriptorBinding("normalMap", "VSG_NORMAL_MAP", MATERIAL_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec3Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32G32B32_SFLOAT}));
        shaderSet->addDescriptorBinding("aoMap", "VSG_LIGHTMAP_MAP", MATERIAL_DESCRIPTOR_SET, 3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::floatArray2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));
        shaderSet->addDescriptorBinding("emissiveMap", "VSG_EMISSIVE_MAP", MATERIAL_DESCRIPTOR_SET, 4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::ubvec4Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R8G8B8A8_UNORM}));
        shaderSet->addDescriptorBinding("material", "", MATERIAL_DESCRIPTOR_SET, 10, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::PhongMaterialValue::create());

        shaderSet->addDescriptorBinding("lightData", "", VIEW_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec4Array::create(64));
        shaderSet->addDescriptorBinding("viewportData", "", VIEW_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec4Value::create(0,0, 1280, 1024));
        shaderSet->addDescriptorBinding("shadowMaps", "", VIEW_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::floatArray3D::create(1, 1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));

        shaderSet->addDescriptorBinding("transmittanceTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("irradianceTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("scatteringTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("singleMieScatteringTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("settings", "", ATMOSHPERE_DESCRIPTOR_SET, 4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("positional", "", POSITIONAL_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->defaultShaderHints = vsg::ShaderCompileSettings::create();

        // additional defines
        shaderSet->optionalDefines = {"ATMOSHPERE_RADIANCE", "ATMOSHPERE_VIEWER_IN_SPACE", "ATMOSHPERE_COMBINED_SCATTERING_TEXTURES", "VSG_GREYSACLE_DIFFUSE_MAP", "VSG_TWO_SIDED_LIGHTING", "VSG_POINT_SPRITE"};
        if(atmosphereBinding->radiance)
            shaderSet->defaultShaderHints->defines.insert("ATMOSHPERE_RADIANCE");

        shaderSet->addPushConstantRange("pc", "", VK_SHADER_STAGE_VERTEX_BIT|VK_SHADER_STAGE_FRAGMENT_BIT, 0, 128);

        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_INSTANCE_POSITIONS", "VSG_DISPLACEMENT_MAP"}, vsg::PositionAndDisplacementMapArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_INSTANCE_POSITIONS"}, vsg::PositionArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_DISPLACEMENT_MAP"}, vsg::DisplacementMapArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_BILLBOARD"}, vsg::BillboardArrayState::create()});

        shaderSet->customDescriptorSetBindings.push_back(vsg::ViewDependentStateBinding::create(VIEW_DESCRIPTOR_SET));
        shaderSet->customDescriptorSetBindings.push_back(positionalBinding);
        shaderSet->customDescriptorSetBindings.push_back(atmosphereBinding);

        phongShaderSet = shaderSet;

        return true;
    }

    bool AtmosphereRuntime::createPBRShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts)
    {
        vsg::info("Local pbr_ShaderSet(",options,")");

        auto vertexModule = options->getRefObject<vsg::ShaderModule>(std::string(pbrVertShader));
        auto fragmentModule = options->getRefObject<vsg::ShaderModule>(std::string(pbrFragShader));

        if (!vertexModule || !fragmentModule)
        {
            vsg::warn("Could not find shaders. PBR support disabled");
            return false;
        }
        else
            vsg::info("Created PBRShaderSet");

        auto vertexShader = vsg::ShaderStage::create(VK_SHADER_STAGE_VERTEX_BIT, "main", vertexModule);
        auto fragmentShader = vsg::ShaderStage::create(VK_SHADER_STAGE_FRAGMENT_BIT, "main", fragmentModule);
        fragmentShader->specializationConstants = constatnts;


        auto shaderSet = vsg::ShaderSet::create(vsg::ShaderStages{vertexShader, fragmentShader});

        shaderSet->addAttributeBinding("vsg_Vertex", "", 0, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_Normal", "", 1, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_TexCoord0", "", 2, VK_FORMAT_R32G32_SFLOAT, vsg::vec2Array::create(1));
        shaderSet->addAttributeBinding("vsg_Color", "", 3, VK_FORMAT_R32G32B32A32_SFLOAT, vsg::vec4Array::create(1));

        shaderSet->addAttributeBinding("vsg_position", "VSG_INSTANCE_POSITIONS", 4, VK_FORMAT_R32G32B32_SFLOAT, vsg::vec3Array::create(1));
        shaderSet->addAttributeBinding("vsg_position_scaleDistance", "VSG_BILLBOARD", 4, VK_FORMAT_R32G32B32A32_SFLOAT, vsg::vec4Array::create(1));

        shaderSet->addDescriptorBinding("displacementMap", "VSG_DISPLACEMENT_MAP", MATERIAL_DESCRIPTOR_SET, 6, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_VERTEX_BIT, vsg::floatArray2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));
        shaderSet->addDescriptorBinding("diffuseMap", "VSG_DIFFUSE_MAP", MATERIAL_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::ubvec4Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R8G8B8A8_UNORM}));
        shaderSet->addDescriptorBinding("mrMap", "VSG_METALLROUGHNESS_MAP", MATERIAL_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec2Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32G32_SFLOAT}));
        shaderSet->addDescriptorBinding("normalMap", "VSG_NORMAL_MAP", MATERIAL_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec3Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32G32B32_SFLOAT}));
        shaderSet->addDescriptorBinding("aoMap", "VSG_LIGHTMAP_MAP", MATERIAL_DESCRIPTOR_SET, 3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::floatArray2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));
        shaderSet->addDescriptorBinding("emissiveMap", "VSG_EMISSIVE_MAP", MATERIAL_DESCRIPTOR_SET, 4, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::ubvec4Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R8G8B8A8_UNORM}));
        shaderSet->addDescriptorBinding("specularMap", "VSG_SPECULAR_MAP", MATERIAL_DESCRIPTOR_SET, 5, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::ubvec4Array2D::create(1, 1, vsg::Data::Properties{VK_FORMAT_R8G8B8A8_UNORM}));
        shaderSet->addDescriptorBinding("material", "", MATERIAL_DESCRIPTOR_SET, 10, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::PbrMaterialValue::create());

        shaderSet->addDescriptorBinding("lightData", "", VIEW_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec4Array::create(64));
        shaderSet->addDescriptorBinding("viewportData", "", VIEW_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, vsg::vec4Value::create(0,0, 1280, 1024));
        shaderSet->addDescriptorBinding("shadowMaps", "", VIEW_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, vsg::floatArray3D::create(1, 1, 1, vsg::Data::Properties{VK_FORMAT_R32_SFLOAT}));

        shaderSet->addDescriptorBinding("transmittanceTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("irradianceTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("scatteringTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("singleMieScatteringTexture", "", ATMOSHPERE_DESCRIPTOR_SET, 3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("settings", "", ATMOSHPERE_DESCRIPTOR_SET, 4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("positional", "", POSITIONAL_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("BRDFTexture", "", PBR_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("specularTexture", "", PBR_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->defaultShaderHints = vsg::ShaderCompileSettings::create();

        // additional defines
        shaderSet->optionalDefines = {"ATMOSHPERE_RADIANCE", "ATMOSHPERE_VIEWER_IN_SPACE", "ATMOSHPERE_COMBINED_SCATTERING_TEXTURES", "VSG_GREYSCALE_DIFFUSE_MAP", "VSG_TWO_SIDED_LIGHTING", "VSG_WORKFLOW_SPECGLOSS"};
        if(atmosphereBinding->radiance)
            shaderSet->defaultShaderHints->defines.insert("ATMOSHPERE_RADIANCE");

        shaderSet->addPushConstantRange("pc", "", VK_SHADER_STAGE_VERTEX_BIT|VK_SHADER_STAGE_FRAGMENT_BIT, 0, 128);

        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_INSTANCE_POSITIONS", "VSG_DISPLACEMENT_MAP"}, vsg::PositionAndDisplacementMapArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_INSTANCE_POSITIONS"}, vsg::PositionArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_DISPLACEMENT_MAP"}, vsg::DisplacementMapArrayState::create()});
        shaderSet->definesArrayStates.push_back(vsg::DefinesArrayState{{"VSG_BILLBOARD"}, vsg::BillboardArrayState::create()});

        shaderSet->customDescriptorSetBindings.push_back(vsg::ViewDependentStateBinding::create(VIEW_DESCRIPTOR_SET));
        shaderSet->customDescriptorSetBindings.push_back(positionalBinding);
        shaderSet->customDescriptorSetBindings.push_back(atmosphereBinding);
        shaderSet->customDescriptorSetBindings.push_back(pbrBinding);

        pbrShaderSet = shaderSet;

        return true;
    }

    bool AtmosphereRuntime::createSkyShaderSet(vsg::ref_ptr<vsg::Options> options, const vsg::ShaderStage::SpecializationConstants &constatnts)
    {
        auto vertexModule = options->getRefObject<vsg::ShaderModule>(std::string(skyVertShader));
        auto fragmentModule = options->getRefObject<vsg::ShaderModule>(std::string(skyFragShader));

        if (!vertexModule || !fragmentModule)
        {
            vsg::error("SkyShaderSet(...) could not find shaders.");
            return false;
        }
        auto vertexShader = vsg::ShaderStage::create(VK_SHADER_STAGE_VERTEX_BIT, "main", vertexModule);
        auto fragmentShader = vsg::ShaderStage::create(VK_SHADER_STAGE_FRAGMENT_BIT, "main", fragmentModule);

        fragmentShader->specializationConstants = constatnts;

        auto shaderSet = vsg::ShaderSet::create(vsg::ShaderStages{vertexShader, fragmentShader});

        shaderSet->addDescriptorBinding("s_Transmittance", "", ATMOSHPERE_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("s_Irradiance", "", ATMOSHPERE_DESCRIPTOR_SET, 1, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("s_Scattering", "", ATMOSHPERE_DESCRIPTOR_SET, 2, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});
        shaderSet->addDescriptorBinding("s_SingleMieScattering", "", ATMOSHPERE_DESCRIPTOR_SET, 3, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("settings", "", ATMOSHPERE_DESCRIPTOR_SET, 4, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->addDescriptorBinding("positional", "", POSITIONAL_DESCRIPTOR_SET, 0, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1, VK_SHADER_STAGE_FRAGMENT_BIT, {});

        shaderSet->defaultShaderHints = vsg::ShaderCompileSettings::create();

        shaderSet->optionalDefines = {"ATMOSHPERE_RADIANCE", "ATMOSHPERE_VIEWER_IN_SPACE"};
        if(atmosphereBinding->radiance)
            shaderSet->defaultShaderHints->defines.insert("ATMOSHPERE_RADIANCE");

        shaderSet->addPushConstantRange("pc", "", VK_SHADER_STAGE_VERTEX_BIT|VK_SHADER_STAGE_FRAGMENT_BIT, 0, 128);

        shaderSet->customDescriptorSetBindings.push_back(atmosphereBinding);
        shaderSet->customDescriptorSetBindings.push_back(inversePositionalBinding);
        if(cloudsBinding)
            shaderSet->customDescriptorSetBindings.push_back(cloudsBinding);

        skyShaderSet = shaderSet;

        return true;
    }

    void AtmosphereRuntime::read(vsg::Input &input)
    {
        //input.read("cubeSize", cubeSize);
        input.read("numViewerThreads", numViewerThreads);

        input.read("lengthUnitInMeters", lengthUnitInMeters);
        /*
            input.readObject("transmittanceData", transmittanceData);
            input.readObject("irradianceData", irradianceData);
            input.readObject("scatteringData", scatteringData);
            input.readObject("singleMieScatteringData", singleMieScatteringData);

            auto sampler = vsg::Sampler::create();
            transmittanceTexture = vsg::ImageInfo::create(sampler, transmittanceData);
            irradianceTexture = vsg::ImageInfo::create(sampler, irradianceData);
            scatteringTexture = vsg::ImageInfo::create(sampler, scatteringData);
            singleMieScatteringTexture = vsg::ImageInfo::create(sampler, singleMieScatteringData);
        */
        //input.readObject("reflectionMapShader", reflectionMapShader);
        //input.readObject("environmentMapShader", environmentMapShader);
        input.readObject("phongShaderSet", phongShaderSet);
        input.readObject("pbrShaderSet", pbrShaderSet);

        input.readObject("ellipsoidModel", ellipsoidModel);

        //auto runtime = input.readValue<RuntimeSettings>("runtimeSettings");

        //runtimeSettings = vsg::Value<RuntimeSettings>::create(runtime);
        //runtimeSettings->properties.dataVariance = vsg::DYNAMIC_DATA;
    }

    void AtmosphereRuntime::write(vsg::Output &output) const
    {
        //output.write("cubeSize", cubeSize);
        output.write("numViewerThreads", numViewerThreads);

        output.write("lengthUnitInMeters", lengthUnitInMeters);
        /*
            output.writeObject("transmittanceData", transmittanceData);
            output.writeObject("irradianceData", irradianceData);
            output.writeObject("scatteringData", scatteringData);
            output.writeObject("singleMieScatteringData", singleMieScatteringData);

            output.writeObject("reflectionMapShader", reflectionMapShader);
            output.writeObject("environmentMapShader", environmentMapShader);
            output.writeObject("phongShaderSet", phongShaderSet);
            output.writeObject("pbrShaderSet", pbrShaderSet);
        */
        output.writeObject("ellipsoidModel", ellipsoidModel);

        //output.write("runtimeSettings", runtimeSettings->value());
    }

    void AtmosphereRuntime::setSunAngle(double radians)
    {
        sunDirection = {0.0, std::sin(radians + vsg::PI), std::cos(radians + vsg::PI)};
    }

    void AtmosphereRuntime::setDate(tm time)
    {
        // Set time(s) relative to the J2000.0 epoch
        //fractionalDay = 367 * y - 7 * (y + (m + 9) \ 12) \ 4 + 275 * m \ 9 + d - 730531.5 + h / 24
        auto fractionalDay = day2000(time);
        // Mean longitude of the Sun
        auto meanLongitudeSunDegrees = 280.4606184 + ((36000.77005361 / 36525) * fractionalDay); //(degrees)
        // Mean anomaly of the Sun
        auto meanAnomalySunDegrees = 357.5277233 +
                                     ((35999.05034 / 36525) * fractionalDay); //(degrees)
        auto meanAnomalySunRadians = meanAnomalySunDegrees * vsg::PI / 180;
        // Ecliptic longitude of the Sun
        auto eclipticLongitudeSunDegrees = meanLongitudeSunDegrees +
                                           (1.914666471 * sin(meanAnomalySunRadians)) +
                                           (0.918994643 * sin(2 * meanAnomalySunRadians)); //(degrees)
        auto eclipticLongitudeSunRadians = eclipticLongitudeSunDegrees * vsg::PI / 180;
        // Obliquity of the ecliptic plane formula mostly derived from:
        // https://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic
        auto epsilonDegrees = // Formula deals with time denominated in centuries
            23.43929 - ((46.8093 / 3600) * fractionalDay / 36525);  //(degrees)
        auto epsilonRadians = epsilonDegrees * vsg::PI / 180;

        vsg::dvec3 direction{cos(eclipticLongitudeSunRadians),
                             cos(epsilonRadians) * sin(eclipticLongitudeSunRadians),
                             sin(epsilonRadians) * sin(eclipticLongitudeSunRadians)};

        sunDirection = direction;
    }


    vsg::ref_ptr<vsg::CommandGraph> AtmosphereRuntime::createCubeMapGraph(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Options> options)
    {
        auto shaderModule = options->getRefObject<vsg::ShaderModule>(std::string(prefilterShader));
        if (!shaderModule)
        {
            vsg::error("Prefilter cubemap(...) could not find shaders.");
            return {};
        }
        shaderModule->hints = vsg::ShaderCompileSettings::create();
        if(atmosphereBinding->radiance)
            shaderModule->hints->defines.insert("ATMOSPHERE_RADIANCE");

        auto size = pbrBinding->prefilteredTexture->imageView->image->extent.width;

        auto shaderStage = vsg::ShaderStage::create(VK_SHADER_STAGE_COMPUTE_BIT, "main", shaderModule);
        shaderStage->specializationConstants.insert_or_assign(1, vsg::intValue::create(numViewerThreads));
        shaderStage->specializationConstants.insert_or_assign(2, vsg::intValue::create(size));

        vsg::DescriptorSetLayoutBindings descriptorBindings{
            {0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1, VK_SHADER_STAGE_COMPUTE_BIT, nullptr}
        };
        auto descriptorSetLayout = vsg::DescriptorSetLayout::create(descriptorBindings);

        vsg::PushConstantRanges pushConstantRanges{
            //{VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(vsg::vec4)}
        };

        auto descriptorImage = vsg::DescriptorImage::create(pbrBinding->prefilteredTexture, 0, 0, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE);

        vsg::DescriptorSetLayouts descriptorSetLayouts{
            vsg::DescriptorSetLayout::create(),
            vsg::DescriptorSetLayout::create(),
            atmosphereBinding->createDescriptorSetLayout(),
            positionalBinding->createDescriptorSetLayout(),
            vsg::DescriptorSetLayout::create(),
            descriptorSetLayout
        };
        auto pipelineLayout = vsg::PipelineLayout::create(descriptorSetLayouts, pushConstantRanges);

        auto descriptorSet = vsg::DescriptorSet::create(descriptorSetLayout, vsg::Descriptors{descriptorImage});
        auto bindDescriptorSet = vsg::BindDescriptorSet::create(VK_PIPELINE_BIND_POINT_COMPUTE, pipelineLayout, descriptorSet);
        auto bindAtmosphere = atmosphereBinding->createStateCommand(pipelineLayout);
        auto bindPositional = positionalBinding->createStateCommand(pipelineLayout);

        // set up the compute pipeline
        auto pipeline = vsg::ComputePipeline::create(pipelineLayout, shaderStage);
        auto bindPipeline = vsg::BindComputePipeline::create(pipeline);
/*
        auto preCopyBarrier = vsg::ImageMemoryBarrier::create();
        preCopyBarrier->srcAccessMask = 0;
        preCopyBarrier->dstAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
        preCopyBarrier->oldLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        preCopyBarrier->newLayout = VK_IMAGE_LAYOUT_GENERAL;
        preCopyBarrier->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        preCopyBarrier->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        preCopyBarrier->image = pbrBinding->prefilteredTexture->imageView->image;
        preCopyBarrier->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        preCopyBarrier->subresourceRange.baseArrayLayer = 0;
        preCopyBarrier->subresourceRange.layerCount = 6;
        preCopyBarrier->subresourceRange.levelCount = 1;
        preCopyBarrier->subresourceRange.baseMipLevel = 0;

        auto preCopyBarrierCmd = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, preCopyBarrier);

        auto postCopyBarrier = vsg::ImageMemoryBarrier::create();
        postCopyBarrier->srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
        postCopyBarrier->dstAccessMask = VK_ACCESS_SHADER_READ_BIT;
        postCopyBarrier->oldLayout = VK_IMAGE_LAYOUT_GENERAL;
        postCopyBarrier->newLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL;
        postCopyBarrier->srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        postCopyBarrier->dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        postCopyBarrier->image = pbrBinding->prefilteredTexture->imageView->image;;
        postCopyBarrier->subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        postCopyBarrier->subresourceRange.baseArrayLayer = 0;
        postCopyBarrier->subresourceRange.layerCount = 6;
        postCopyBarrier->subresourceRange.levelCount = 1;
        postCopyBarrier->subresourceRange.baseMipLevel = 0;

        auto postCopyBarrierCmd = vsg::PipelineBarrier::create(VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, 0, postCopyBarrier);
*/
        int computeQueueFamily = window->getOrCreatePhysicalDevice()->getQueueFamily(VK_QUEUE_COMPUTE_BIT);
        auto compute_commandGraph = vsg::CommandGraph::create(window->getOrCreateDevice(), computeQueueFamily);

        //compute_commandGraph->addChild(preCopyBarrierCmd);
        compute_commandGraph->addChild(bindPipeline);
        compute_commandGraph->addChild(bindAtmosphere);
        compute_commandGraph->addChild(bindPositional);
        compute_commandGraph->addChild(bindDescriptorSet);
        auto workgroups = uint32_t(ceil(float(size) / float(numViewerThreads)));
        compute_commandGraph->addChild(vsg::Dispatch::create(workgroups, workgroups, 6));
        //compute_commandGraph->addChild(postCopyBarrierCmd);

        return compute_commandGraph;
    }

    vsg::ref_ptr<vsg::Node> AtmosphereRuntime::createSky(bool viewerInSpace)
    {
        auto drawCommands = vsg::Commands::create();

        auto indices = vsg::ushortArray::create({0, 2, 1, 1, 2, 3});

        drawCommands->addChild(vsg::BindIndexBuffer::create(indices));
        drawCommands->addChild(vsg::DrawIndexed::create(indices->size(), 1, 0, 0, 0));

        auto graphicsPipelineConfig = vsg::GraphicsPipelineConfigurator::create(skyShaderSet);

        auto& defines = graphicsPipelineConfig->shaderHints->defines;

        if(viewerInSpace)
            defines.insert("ATMOSHPERE_VIEWER_IN_SPACE");

        struct SetPipelineStates : public vsg::Visitor
        {
            void apply(Object& object) override
            {
                object.traverse(*this);
            }

            void apply(vsg::DepthStencilState& dss) override
            {
                dss.depthCompareOp = VK_COMPARE_OP_GREATER_OR_EQUAL;
            }
        };
        SetPipelineStates sps;

        graphicsPipelineConfig->accept(sps);

        graphicsPipelineConfig->init();

        vsg::StateCommands stateCommands;
        if (graphicsPipelineConfig->copyTo(stateCommands))
        {
            // create StateGroup as the root of the scene/command graph to hold the GraphicsPipeline, and binding of Descriptors to decorate the whole graph
            auto stateGroup = vsg::StateGroup::create();
            stateGroup->stateCommands.swap(stateCommands);
            stateGroup->prototypeArrayState = graphicsPipelineConfig->getSuitableArrayState();
            stateGroup->addChild(drawCommands);
            return stateGroup;
        }
        return {};
    }

    vsg::ref_ptr<vsg::View> AtmosphereRuntime::createSkyView(vsg::ref_ptr<vsg::Window> window, vsg::ref_ptr<vsg::Camera> camera)
    {
        // create the sky camera
        auto inversePerojection = InverseProjection::create(camera->projectionMatrix);
        auto inverseView = InverseView::create(camera->viewMatrix);
        auto skyCamera = vsg::Camera::create(inversePerojection, inverseView, vsg::ViewportState::create(window->extent2D()));

        return vsg::View::create(skyCamera);
    }
}
