cmake_minimum_required(VERSION 3.20)
project(PostTopOptVolumeSurface LANGUAGES CXX)
find_package(ParaView REQUIRED COMPONENTS VTKExtensionsAMR)

# VTK operations and the application pipeline are the only two source files.
add_library(PostTopOptVolumeSurface STATIC VolumeSurface.cpp VtkopContract.cpp)
target_compile_features(PostTopOptVolumeSurface PUBLIC cxx_std_17)
target_include_directories(PostTopOptVolumeSurface PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(PostTopOptVolumeSurface PUBLIC
    ParaView::VTKExtensionsAMR VTK::ParallelCore
    VTK::CommonCore VTK::CommonDataModel VTK::CommonExecutionModel
    VTK::FiltersCore VTK::FiltersGeneral VTK::FiltersExtraction)
vtk_module_autoinit(TARGETS PostTopOptVolumeSurface MODULES
    ParaView::VTKExtensionsAMR VTK::ParallelCore
    VTK::FiltersCore VTK::FiltersGeneral VTK::FiltersExtraction)
