set(SOURCES
        vsgatmosphere.cpp)

add_executable(vsgatmosphere ${SOURCES})

target_link_libraries(vsgatmosphere vsg::vsg vsgAtmosphere::vsgAtmosphere)

if (vsgXchange_FOUND)
    target_compile_definitions(vsgatmosphere PRIVATE vsgXchange_FOUND)
    target_link_libraries(vsgatmosphere vsgXchange::vsgXchange)
endif()

install(TARGETS vsgatmosphere RUNTIME DESTINATION bin)
