# Custom configuration for benchmarks
add_library(BenchConfiguration INTERFACE)
target_link_libraries(BenchConfiguration INTERFACE BaseConfiguration)
if(BUILD_TYPE STREQUAL "Release")
    # Add custom optimization flags here:
    target_compile_options(BenchConfiguration INTERFACE -O3)
else()
    target_link_libraries(BenchConfiguration INTERFACE BaseOptConfiguration)
endif()

function(add_perf_test SOURCE)
    get_filename_component(NAME ${SOURCE} NAME_WE)
    set(NAME "${NAME}_perf")
    add_executable(${NAME} ${SOURCE})

    target_include_directories(${NAME}
        PUBLIC
        ${PROJECT_SOURCE_DIR}/lib/include
        ${PROJECT_SOURCE_DIR}/perf/include
    )

    target_link_libraries(${NAME}
    PUBLIC
        base_naive
    PRIVATE
        BenchConfiguration
    )
endfunction()

file(GLOB SOURCES "*/*.cpp")

foreach(SOURCE ${SOURCES})
    add_perf_test(${SOURCE})
endforeach()
