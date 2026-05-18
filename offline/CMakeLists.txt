cmake_minimum_required(VERSION 3.15)
project(BOProject LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)
set(ENABLE_TESTS OFF CACHE BOOL "Disable tests for dependencies" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "Disable examples for dependencies" FORCE)

FetchContent_Declare(
    kfr
    GIT_REPOSITORY https://github.com/kfrlib/kfr.git
    GIT_TAG 7.0.1
    GIT_SHALLOW TRUE
)
FetchContent_Declare(
    json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_Declare(
    limbo
    GIT_REPOSITORY https://github.com/resibots/limbo.git
    GIT_TAG v2.1.0
    GIT_SHALLOW TRUE
)
FetchContent_Declare(
    Neun
    GIT_REPOSITORY https://github.com/gnb-uam/neun.git
    GIT_TAG 437c5b95898a0542f2c96a5344856d8dac12ad52
)
FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
    GIT_TAG yaml-cpp-0.9.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(kfr json limbo Neun yaml-cpp)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem thread)
find_package(TBB REQUIRED)
find_package(NLopt REQUIRED)

add_executable(consts_calculator consts_calculator.cpp scaling.cpp)
target_compile_options(consts_calculator PRIVATE -march=native)
target_include_directories(consts_calculator PRIVATE
    ${neun_SOURCE_DIR}
    ${neun_SOURCE_DIR}/include
    ${neun_SOURCE_DIR}/archetypes
    ${neun_SOURCE_DIR}/concepts
    ${neun_SOURCE_DIR}/integrators
    ${neun_SOURCE_DIR}/models
    ${neun_SOURCE_DIR}/wrappers
)
target_link_libraries(consts_calculator PRIVATE kfr kfr_dsp)

add_executable(bo main.cpp scaling.cpp)
target_compile_options(bo PRIVATE -march=native)
target_compile_definitions(bo PRIVATE USE_NLOPT)
target_include_directories(bo PRIVATE
    ${limbo_SOURCE_DIR}/src
    ${neun_SOURCE_DIR}
    ${neun_SOURCE_DIR}/include
    ${neun_SOURCE_DIR}/archetypes
    ${neun_SOURCE_DIR}/concepts
    ${neun_SOURCE_DIR}/integrators
    ${neun_SOURCE_DIR}/models
    ${neun_SOURCE_DIR}/wrappers
)
target_link_libraries(bo PRIVATE
    kfr
    kfr_dsp
    Eigen3::Eigen
    Boost::system
    Boost::filesystem
    Boost::thread
    TBB::tbb
    nlopt
    yaml-cpp
    nlohmann_json::nlohmann_json
)