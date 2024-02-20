file(GLOB SOURCES "*/*.cpp")

add_executable(most_common_string_test ${SOURCES})

target_include_directories(most_common_string_test
    PUBLIC
    ${PROJECT_SOURCE_DIR}/test/include
    ${PROJECT_SOURCE_DIR}/3rd-party/gtest/googletest/include
    ${PROJECT_SOURCE_DIR}/lib/include
)

target_link_libraries(most_common_string_test
PUBLIC
    GTest::gtest
    GTest::gtest_main
    base_naive
PRIVATE
    CommonConfiguration
)


