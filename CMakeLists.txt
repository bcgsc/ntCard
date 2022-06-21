project(ntCard)
cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

include_directories(include)
include_directories(vendor/argparse/include)
include_directories(vendor/btllib/include)
include_directories(vendor/btllib/subprojects/sdsl-lite/include)
include_directories(vendor/btllib/subprojects/cpptoml/include)

file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a
        COMMAND ${PROJECT_SOURCE_DIR}/vendor/btllib/compile --prefix ${PROJECT_BINARY_DIR}/btllib
        USES_TERMINAL
)
set_property(SOURCE main.cpp APPEND PROPERTY OBJECT_DEPENDS ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)

add_executable(ntCard ${SRC_FILES} main.cpp)
target_link_libraries(ntCard ${PROJECT_BINARY_DIR}/btllib/lib/libbtllib.a)