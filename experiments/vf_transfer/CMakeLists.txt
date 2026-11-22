cmake_minimum_required(VERSION 3.20)
project(PostTopOptVFTransferExperiments LANGUAGES CXX)

# Use the canonical integration component; no algorithm sources live here.
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../integration/vtkop"
                 "${CMAKE_CURRENT_BINARY_DIR}/vtkop")

# Standalone end-to-end shape transfer and AMRDualContour OBJ export.

include(CTest)
if(BUILD_TESTING)
    add_executable(VFTransferShapes ShapesTest.cpp)
    target_link_libraries(VFTransferShapes PRIVATE PostTopOptVolumeSurface ParaView::VTKExtensionsAMR)
    vtk_module_autoinit(TARGETS VFTransferShapes MODULES ParaView::VTKExtensionsAMR VTK::FiltersCore)
    add_test(NAME VFTransferShapes COMMAND VFTransferShapes "${CMAKE_CURRENT_BINARY_DIR}/shape-results")
endif()
# Explicit build target: builds and runs all shapes; never part of ALL.
if(BUILD_TESTING)
    add_custom_target(vf-transfer-shapes
        COMMAND "${CMAKE_COMMAND}" -DTEST_EXE=$<TARGET_FILE:VFTransferShapes> -DOUTPUT_ROOT=${CMAKE_CURRENT_BINARY_DIR} -P "${CMAKE_CURRENT_SOURCE_DIR}/RunShapes.cmake"
        DEPENDS VFTransferShapes
        USES_TERMINAL VERBATIM)
endif()


# Example host application: FE arrays -> public integration API -> Surface -> OBJ.
add_executable(IntegrationExample IntegrationExample.cpp)
target_link_libraries(IntegrationExample PRIVATE PostTopOptVolumeSurface)
