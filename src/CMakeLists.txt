SET(HEADER_PATH ${VSGATMOSPHERE_SOURCE_DIR}/include/vsgAtmosphere)

set(SOURCES
        Atmosphere.cpp
        AtmosphereLighting.cpp
        AtmosphereBinding.cpp
        AtmosphereRuntime.cpp
        AtmosphereImage.cpp
        AtmosphereTrackball.cpp
        Clouds.cpp
)
set(HEADERS
        ${HEADER_PATH}/AtmoshpereConstatnts.h
        ${HEADER_PATH}/AtmosphereTools.h
        ${HEADER_PATH}/Atmosphere.h
        ${HEADER_PATH}/AtmosphereLighting.h
        ${HEADER_PATH}/AtmosphereBinding.h
        ${HEADER_PATH}/AtmosphereRuntime.h
        ${HEADER_PATH}/InverseMatrices.h
        ${HEADER_PATH}/AtmosphereImage.h
        ${HEADER_PATH}/AtmosphereTrackball.h
        ${HEADER_PATH}/Clouds.h
)

add_library(vsgAtmosphere ${HEADERS} ${SOURCES})

target_include_directories(vsgAtmosphere
    PUBLIC
        $<BUILD_INTERFACE:${HEADER_PATH}>
        $<BUILD_INTERFACE:${HEADER_PATH}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(vsgAtmosphere PUBLIC vsg::vsg vsgXchange::vsgXchange)

# add definitions to enable building vsgAtmosphere as part of submodule
add_library(vsgAtmosphere::vsgAtmosphere ALIAS vsgAtmosphere)
set(vsgAtmosphere_FOUND TRUE CACHE INTERNAL "vsgAtmosphere found.")
set(CMAKE_DISABLE_FIND_PACKAGE_vsgAtmosphere TRUE CACHE INTERNAL "Disable find_package(vsgAtmosphere) as it's not necessary.")


set_property(TARGET vsgAtmosphere PROPERTY VERSION ${VSGXCHANGE_VERSION_MAJOR}.${VSGXCHANGE_VERSION_MINOR}.${VSGXCHANGE_VERSION_PATCH})
set_property(TARGET vsgAtmosphere PROPERTY SOVERSION ${VSGXCHANGE_SOVERSION})
set_property(TARGET vsgAtmosphere PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET vsgAtmosphere PROPERTY CXX_STANDARD 17)

target_compile_definitions(vsgAtmosphere PRIVATE ${EXTRA_DEFINES})

install(TARGETS vsgAtmosphere ${INSTALL_TARGETS_DEFAULT_FLAGS})

install(DIRECTORY include/vsgXchange DESTINATION include)
if (NOT(${VSGATMOSPHERE_BINARY_DIR} STREQUAL ${VSGATMOSPHERE_SOURCE_DIR}))
    install(DIRECTORY ${VSGATMOSPHERE_BINARY_DIR}/include/vsgXchange DESTINATION include)
endif()

string(REPLACE ";" "\n" FIND_DEPENDENCY_OUT "${FIND_DEPENDENCY}")
