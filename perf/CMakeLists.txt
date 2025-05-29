# Custom configuration for benchmarks
add_library(BenchConfiguration INTERFACE)
set(CMAKE_EXE_LINKER_FLAGS "-static ${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_FIND_LIBRARY_SUFFIXES .a)  # Искать только статические библиотеки
target_link_libraries(BenchConfiguration INTERFACE BaseConfiguration)
if(BUILD_TYPE STREQUAL "Release")
    # Add custom optimization flags here:
    target_compile_options(BenchConfiguration INTERFACE -Ofast)
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
            Kokkos::kokkos
        PRIVATE
            BenchConfiguration
            ${CMAKE_DL_LIBS}
            pthread
    )
    
    # Для статической линковки
    
        target_link_options(${NAME} PRIVATE -static)
    
endfunction()

file(GLOB SOURCES "*/*.cpp")

foreach(SOURCE ${SOURCES})
    add_perf_test(${SOURCE})
endforeach()
