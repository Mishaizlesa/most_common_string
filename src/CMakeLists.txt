set(target ${MP2_TESTS})

file(GLOB hdrs "*.h*")
file(GLOB srcs "*.cpp")


file(GLOB samples_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(sample_filename ${samples_list})
  # Get file name without extension
  get_filename_component(sample ${sample_filename} NAME_WE)

  # Add and configure executable file to be produced
  add_executable(${sample} ${sample_filename})
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
      target_link_libraries(${sample} OpenMP::OpenMP_CXX)
  endif()
endforeach()