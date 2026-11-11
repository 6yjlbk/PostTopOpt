
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
    GIT_TAG 3.4.1 ## Версия Eigen 3.4.1
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(Eigen)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(PostTopOpt STATIC src/DualContouring.cpp src/MarchingCubes.cpp src/ObjWriter.cpp src/CgalExtraction.cpp)
target_include_directories(PostTopOpt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(PostTopOpt PUBLIC Eigen3::Eigen)
add_executable(PostTopOptDemo examples/sphere.cpp)
target_link_libraries(PostTopOptDemo PRIVATE PostTopOpt)
target_compile_definitions(PostTopOptDemo PRIVATE
    POSTTOPOPT_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results")
include(CTest)
if(BUILD_TESTING)
    add_executable(PostTopOptTest tests/dual_contouring.cpp)
    target_link_libraries(PostTopOptTest PRIVATE PostTopOpt)
    add_test(NAME DualContouringTest COMMAND PostTopOptTest)
endif()

# Analytic geometry examples; keep the original sphere executable unchanged.
foreach(example IN ITEMS Torus SphereUnion DrilledBox Gyroid)
    if(example STREQUAL "Torus")
        set(example_source torus)
    elseif(example STREQUAL "SphereUnion")
        set(example_source sphere_union)
    elseif(example STREQUAL "DrilledBox")
        set(example_source drilled_box)
    else()
        set(example_source gyroid)
    endif()
    add_executable(PostTopOpt${example} examples/${example_source}.cpp)
    target_link_libraries(PostTopOpt${example} PRIVATE PostTopOpt)
target_compile_definitions(PostTopOpt${example} PRIVATE
    POSTTOPOPT_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results")
    if(BUILD_TESTING)
        add_test(NAME Example${example}
                 COMMAND PostTopOpt${example} "${CMAKE_CURRENT_SOURCE_DIR}/results/${example_source}.obj")
    endif()
endforeach()

if(BUILD_TESTING)
    add_executable(PostTopOptMethodsTest tests/extraction_methods.cpp)
    target_link_libraries(PostTopOptMethodsTest PRIVATE PostTopOpt)
    add_test(NAME ExtractionMethodsTest COMMAND PostTopOptMethodsTest)
    foreach(shape IN ITEMS SphereUnion DrilledBox Gyroid)
        foreach(method IN ITEMS dc emilk mc)
            add_test(NAME Compare${shape}_${method}
                COMMAND PostTopOpt${shape} "${CMAKE_CURRENT_SOURCE_DIR}/results/${shape}_${method}.obj" --method ${method})
        endforeach()
    endforeach()
    add_test(NAME InvalidMethodCLI COMMAND PostTopOptSphereUnion --method invalid)
    set_tests_properties(InvalidMethodCLI PROPERTIES WILL_FAIL TRUE)
endif()

option(POSTTOPOPT_WITH_CGAL "Enable the optional CGAL DC backend" OFF)
if(POSTTOPOPT_WITH_CGAL)
    find_package(CGAL 6.1 REQUIRED)
    target_link_libraries(PostTopOpt PRIVATE CGAL::CGAL)
    target_compile_definitions(PostTopOpt PRIVATE POSTTOPOPT_WITH_CGAL=1)
endif()

if(BUILD_TESTING)
    add_executable(PostTopOptCgalTest tests/cgal_extraction.cpp)
    target_link_libraries(PostTopOptCgalTest PRIVATE PostTopOpt)
    add_test(NAME CgalBackendTest COMMAND PostTopOptCgalTest)
    if(POSTTOPOPT_WITH_CGAL)
        foreach(shape IN ITEMS SphereUnion DrilledBox Gyroid)
            add_test(NAME CompareCGAL_${shape}
                COMMAND PostTopOpt${shape} "${CMAKE_CURRENT_SOURCE_DIR}/results/${shape}_cgal-dc.obj" --method cgal-dc)
        endforeach()
    endif()
endif()

# Manual four-run gradient experiment; intentionally not registered with CTest.
add_executable(PostTopOptGradientComparison examples/gradient_comparison.cpp)
target_link_libraries(PostTopOptGradientComparison PRIVATE PostTopOpt)
target_compile_definitions(PostTopOptGradientComparison PRIVATE
    POSTTOPOPT_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results")

# Cone examples share the standard CLI; select emilk with --method emilk.
foreach(mode IN ITEMS Numeric Analytic)
    add_executable(PostTopOptCone${mode} examples/cone.cpp)
    target_link_libraries(PostTopOptCone${mode} PRIVATE PostTopOpt)
    target_compile_definitions(PostTopOptCone${mode} PRIVATE
        POSTTOPOPT_RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/results")
endforeach()
target_compile_definitions(PostTopOptConeAnalytic PRIVATE POSTTOPOPT_CONE_ANALYTIC=1)
