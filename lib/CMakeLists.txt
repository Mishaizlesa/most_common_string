set(LIBRARY_NAME "base_naive")

# Создаем основную библиотеку (C++)
add_library(${LIBRARY_NAME}
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/base_naive.cpp"
)
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/include)
target_link_libraries(${LIBRARY_NAME} PRIVATE CommonConfiguration)

# Новая библиотека для Coarray Fortran
set(FORTRAN_LIBRARY_NAME "algo_coarray")
add_library(${FORTRAN_LIBRARY_NAME} STATIC
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/rabin_karp_SWAR.f90" # Добавляем Fortran-файл
)
target_compile_options(${FORTRAN_LIBRARY_NAME} PRIVATE -coarray=shared)
target_include_directories(${FORTRAN_LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/include)
target_link_libraries(${FORTRAN_LIBRARY_NAME} PRIVATE -coarray)

# Модифицированная функция для обработки C++ и Fortran
function(add_obj_lib OBJ_NAME CONFIGURATION)
    # Проверяем расширение файла (.cpp или .f90)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/${OBJ_NAME}.cpp")
        add_library(${OBJ_NAME} OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/${OBJ_NAME}.cpp")
        target_include_directories(${OBJ_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/include)
        target_link_libraries(${OBJ_NAME} PRIVATE ${CONFIGURATION} -fsycl)
        target_link_libraries(${LIBRARY_NAME} PRIVATE ${OBJ_NAME})
    elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/${OBJ_NAME}.f90")
        add_library(${OBJ_NAME} OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/lib/src/${OBJ_NAME}.f90")
        target_include_directories(${OBJ_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/include)
        target_compile_options(${OBJ_NAME} PRIVATE -coarray=shared)
        target_link_libraries(${OBJ_NAME} PRIVATE ${CONFIGURATION} -coarray)
        target_link_libraries(${FORTRAN_LIBRARY_NAME} PRIVATE ${OBJ_NAME})
    else()
        message(FATAL_ERROR "Source file for ${OBJ_NAME} not found")
    endif()
endfunction()

# Кастомная конфигурация
add_library(CustomConfiguration INTERFACE)
target_link_libraries(CustomConfiguration INTERFACE BaseConfiguration)
if(BUILD_TYPE STREQUAL "Release")
    target_compile_options(CustomConfiguration INTERFACE -Ofast)
else()
    target_link_libraries(CustomConfiguration INTERFACE BaseOptConfiguration)
endif()

# Объектные библиотеки
add_obj_lib("naive" CommonConfiguration)
add_obj_lib("hash3" CommonConfiguration)
add_obj_lib("rabin_karp_rolling_hash" CommonConfiguration)
add_obj_lib("rabin_karp_SWAR" CommonConfiguration) # Может быть Fortran или C++