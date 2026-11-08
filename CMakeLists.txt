
cmake_minimum_required(VERSION 3.20)

project(PostTopOpt
    VERSION 0.1.0
    LANGUAGES CXX
)

include(FetchContent)

set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
set(EIGEN_INSTALL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.1
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(Eigen)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(PostTopOpt STATIC src/DualContouring.cpp src/ObjWriter.cpp)
target_include_directories(PostTopOpt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(PostTopOpt PUBLIC Eigen3::Eigen)
add_executable(PostTopOptDemo examples/sphere.cpp)
target_link_libraries(PostTopOptDemo PRIVATE PostTopOpt)
include(CTest)
if(BUILD_TESTING)
    add_executable(PostTopOptTest tests/dual_contouring.cpp)
    target_link_libraries(PostTopOptTest PRIVATE PostTopOpt)
    add_test(NAME DualContouringTest COMMAND PostTopOptTest)
endif()
