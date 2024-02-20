set(LIBRARY_NAME "base_naive")

# Create main library
add_library(${LIBRARY_NAME}
STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src/base_naive.cpp"
)

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(${LIBRARY_NAME} PRIVATE CommonConfiguration)

function(add_obj_lib OBJ_NAME CONFIGURATION)
    add_library(${OBJ_NAME} OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/src/${OBJ_NAME}.cpp")
    target_include_directories(${OBJ_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(${OBJ_NAME} PRIVATE ${CONFIGURATION})
    target_link_libraries(${LIBRARY_NAME} PRIVATE ${OBJ_NAME})
endfunction()

# Custom configuration for separate files example
add_library(CustomConfiguration INTERFACE)
target_link_libraries(CustomConfiguration INTERFACE BaseConfiguration)
if(BUILD_TYPE STREQUAL "Release")
    # Add custom optimization flags here:
    target_compile_options(CustomConfiguration INTERFACE -O3)
else()
    target_link_libraries(CustomConfiguration INTERFACE BaseOptConfiguration)
endif()

#
# Object libraries with separate compiler keys
#
add_obj_lib("naive" CommonConfiguration)
add_obj_lib("hash3" CommonConfiguration)
add_obj_lib("rabin_karp" CommonConfiguration)