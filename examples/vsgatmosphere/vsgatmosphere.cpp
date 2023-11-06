#include "AtmosphereLighting.h"
#include <vsg/all.h>

#ifdef vsgXchange_FOUND
#    include <vsgXchange/all.h>
#endif

#include <algorithm>
#include <chrono>
#include <iostream>
#include <thread>

#include "Atmosphere.h"
#include "InverseMatrices.h"

int main(int argc, char** argv)
{
    try
    {
        // set up defaults and read command line arguments to override them
        vsg::CommandLine arguments(&argc, argv);

        // set up vsg::Options to pass in filepaths and ReaderWriter's and other IO related options to use when reading and writing files.
        auto options = vsg::Options::create();
        options->sharedObjects = vsg::SharedObjects::create();
        options->fileCache = vsg::getEnv("VSG_FILE_CACHE");
        options->paths = vsg::getEnvPaths("VSG_FILE_PATH");

#ifdef vsgXchange_all
        // add vsgXchange's support for reading and writing 3rd party file formats
        options->add(vsgXchange::all::create());
#endif

        arguments.read(options);

        options->paths.emplace_back(options->paths.back() / "shaders" / "scattering");

        auto windowTraits = vsg::WindowTraits::create();
        windowTraits->windowTitle = "vsgatmosphere";
        windowTraits->debugLayer = arguments.read({"--debug", "-d"});
        windowTraits->apiDumpLayer = arguments.read({"--api", "-a"});
        windowTraits->synchronizationLayer = arguments.read("--sync");
        if (int mt = 0; arguments.read({"--memory-tracking", "--mt"}, mt)) vsg::Allocator::instance()->setMemoryTracking(mt);
        if (arguments.read("--double-buffer")) windowTraits->swapchainPreferences.imageCount = 2;
        if (arguments.read("--triple-buffer")) windowTraits->swapchainPreferences.imageCount = 3; // default
        if (arguments.read("--IMMEDIATE")) windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
        if (arguments.read("--FIFO")) windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_FIFO_KHR;
        if (arguments.read("--FIFO_RELAXED")) windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_FIFO_RELAXED_KHR;
        if (arguments.read("--MAILBOX")) windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_MAILBOX_KHR;
        if (arguments.read({"-t", "--test"}))
        {
            windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
            windowTraits->fullscreen = true;
        }
        if (arguments.read({"--st", "--small-test"}))
        {
            windowTraits->swapchainPreferences.presentMode = VK_PRESENT_MODE_IMMEDIATE_KHR;
            windowTraits->width = 192, windowTraits->height = 108;
            windowTraits->decoration = false;
        }
        if (arguments.read({"--fullscreen", "--fs"})) windowTraits->fullscreen = true;
        if (arguments.read({"--window", "-w"}, windowTraits->width, windowTraits->height)) { windowTraits->fullscreen = false; }
        if (arguments.read({"--no-frame", "--nf"})) windowTraits->decoration = false;
        if (arguments.read("--or")) windowTraits->overrideRedirect = true;
        if (arguments.read("--d32")) windowTraits->depthFormat = VK_FORMAT_D32_SFLOAT;
        arguments.read("--screen", windowTraits->screenNum);
        arguments.read("--display", windowTraits->display);
        arguments.read("--samples", windowTraits->samples);
        auto numFrames = arguments.value(-1, "-f");
        auto pathFilename = arguments.value(std::string(), "-p");
        //auto horizonMountainHeight = arguments.value(0.0, "--hmh");

        auto sunAngle = vsg::radians(arguments.value(0.0f, "--angle"));
        auto skyExposure = arguments.value(1.0f, "--sky");
        auto modelExposure = arguments.value(1.0f, "--model");

        if (arguments.read("--rgb")) options->mapRGBtoRGBAHint = false;

        bool generateDebug = arguments.read({"--shader-debug-info", "--sdi"});
        if (generateDebug)
            windowTraits->deviceExtensionNames.push_back(VK_KHR_SHADER_NON_SEMANTIC_INFO_EXTENSION_NAME);

        if (int log_level = 0; arguments.read("--log-level", log_level)) vsg::Logger::instance()->level = vsg::Logger::Level(log_level);

        if (arguments.errors()) return arguments.writeErrorMessages(std::cerr);

        // create the viewer and assign window(s) to it
        auto viewer = vsg::Viewer::create();
        auto window = vsg::Window::create(windowTraits);
        if (!window)
        {
            std::cout << "Could not create windows." << std::endl;
            return 1;
        }
        auto vsg_scene = vsg::Group::create();

        auto ellipsoidModel = vsg::EllipsoidModel::create(vsg::WGS_84_RADIUS_EQUATOR, vsg::WGS_84_RADIUS_EQUATOR);

        viewer->addWindow(window);

        double nearFarRatio = 0.0001;
        auto time = vsg::floatValue::create();
        time->properties.dataVariance = vsg::DYNAMIC_DATA;

        auto modelView = vsg::LookAt::create();
        modelView->center = ellipsoidModel->convertLatLongAltitudeToECEF({51.50151088842245, -0.14181489107549874, 0.0});
        modelView->up = vsg::dvec3(0.0, 1.0, 0.0);
        modelView->eye = vsg::dvec3(40000000.0, 0.0, 0.0);
        //auto perspective = vsg::EllipsoidPerspective::create(modelView, ellipsoidModel, 30.0, static_cast<double>(window->extent2D().width) / static_cast<double>(window->extent2D().height), nearFarRatio, horizonMountainHeight);
        auto perspective = vsg::Perspective::create(60.0, static_cast<double>(window->extent2D().width) / static_cast<double>(window->extent2D().height), nearFarRatio, ellipsoidModel->radiusEquator() * 10.0);

        auto camera = vsg::Camera::create(perspective, modelView, vsg::ViewportState::create(window->extent2D()));

        // create the sky camera
        auto inversePerojection = atmosphere::InverseProjection::create(camera->projectionMatrix);
        auto inverseView = atmosphere::InverseView::create(camera->viewMatrix);
        auto skyCamera = vsg::Camera::create(inversePerojection, inverseView, vsg::ViewportState::create(window->extent2D()));

        vsg::RegisterWithObjectFactoryProxy<atmosphere::AtmosphereModelSettings>();
        vsg::RegisterWithObjectFactoryProxy<atmosphere::AtmosphereData>();

        auto settings = atmosphere::AtmosphereModelSettings::create(ellipsoidModel);

        auto atmosphereGenerator = atmosphere::createAtmosphereGenerator(window, settings, options);

        auto mainViewDependent = atmosphere::AtmosphereLighting::create(modelView);
        atmosphereGenerator->viewDescriptorSetLayout = mainViewDependent->descriptorSetLayout;
        mainViewDependent->exposure = modelExposure;
        auto skyViewDependent = atmosphere::AtmosphereLighting::create(inverseView);
        skyViewDependent->exposure = skyExposure;
        skyViewDependent->transform = false;

        auto atmosphere = atmosphereGenerator->loadData();

        //auto clouds = atmosphere::loadClouds("textures/scattering", options);

        mainViewDependent->assignData(atmosphere);
        skyViewDependent->assignData(atmosphere);

        atmosphere->setSunAngle(sunAngle);

        options->shaderSets["phong"] = atmosphere->phongShaderSet;

        auto builder = vsg::Builder::create();
        builder->options = options;
        auto transfrom = vsg::AbsoluteTransform::create(vsg::translate(0.0, 0.0, -3.0));
        transfrom->addChild(builder->createSphere());
        vsg_scene->addChild(transfrom);

        auto databaseSettings = vsg::createOpenStreetMapSettings(options);
        databaseSettings->lighting = true;
        databaseSettings->ellipsoidModel = ellipsoidModel;

        auto earth = vsg::TileDatabase::create();
        earth->settings = databaseSettings;
        earth->readDatabase(options);

        vsg_scene->addChild(earth);

        // add close handler to respond the close window button and pressing escape
        viewer->addEventHandler(vsg::CloseHandler::create(viewer));

        auto trackball = vsg::Trackball::create(camera, ellipsoidModel);

        trackball->addKeyViewpoint(vsg::KeySymbol('1'), 51.50151088842245, -0.14181489107549874, 10.0, 2.0); // Grenwish Observatory

        viewer->addEventHandler(trackball);

        //auto compute_commandGraph = clouds->createCloudMapGraph(window, mainViewDependent, time);

        // set up the render graph
        //auto renderGraph = vsg::createRenderGraphForView(window, camera, vsg_scene, VK_SUBPASS_CONTENTS_INLINE, false);

        auto mainView = vsg::View::create(camera);
        mainView->addChild(vsg_scene);
        mainView->viewDependentState = mainViewDependent;

        // set up the render graph
        auto renderGraph = vsg::RenderGraph::create(window, mainView);
        renderGraph->contents = VK_SUBPASS_CONTENTS_INLINE;

        auto skyView = vsg::View::create(skyCamera, atmosphere->createSky(skyViewDependent->descriptorSetLayout));
        skyView->viewDependentState = skyViewDependent;

        renderGraph->addChild(skyView);
        renderGraph->setClearValues({{0.0f, 0.0f, 0.0f, 1.0f}});

        auto grahics_commandGraph = vsg::CommandGraph::create(window, renderGraph);
        viewer->assignRecordAndSubmitTaskAndPresentation({grahics_commandGraph});

        viewer->compile();

        // rendering main loop
        while (viewer->advanceToNextFrame() && (numFrames < 0 || (numFrames--) > 0))
        {
            // pass any events into EventHandlers assigned to the Viewer
            viewer->handleEvents();

            viewer->update();

            viewer->recordAndSubmit();

            viewer->present();
        }
    }
    catch (const vsg::Exception& ve)
    {
        for (int i = 0; i < argc; ++i) std::cerr << argv[i] << " ";
        std::cerr << "\n[Exception] - " << ve.message << " result = " << ve.result << std::endl;
        return 1;
    }

    // clean up done automatically thanks to ref_ptr<>
    return 0;
}
