set(SKYSPACES_EULER3D_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SKYSPACES_EULER3D_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SKYSPACES_EULER3D_DEMO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/demo")

add_library(skyspaces_euler3d_reference STATIC
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/AdvectionIntegrator3D.h"
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/AdvectionScheme3D.h"
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/EulerFluid3D.h"
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/Interpolation3D.h"
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/PressureSolver3D.h"
    "${SKYSPACES_EULER3D_INCLUDE_DIR}/Types3D.h"
    "${SKYSPACES_EULER3D_SRC_DIR}/AdvectionIntegrator3D.cpp"
    "${SKYSPACES_EULER3D_SRC_DIR}/AdvectionScheme3D.cpp"
    "${SKYSPACES_EULER3D_SRC_DIR}/EulerFluid3D.cpp"
    "${SKYSPACES_EULER3D_SRC_DIR}/Interpolation3D.cpp"
    "${SKYSPACES_EULER3D_SRC_DIR}/PressureSolver3D.cpp"
)

target_include_directories(skyspaces_euler3d_reference
    PUBLIC
        "${SKYSPACES_EULER3D_INCLUDE_DIR}"
        "${SKYSPACES_INCLUDE_DIR}"
)

target_link_libraries(skyspaces_euler3d_reference
    PUBLIC
        Eigen3::Eigen
)

target_compile_features(skyspaces_euler3d_reference PUBLIC cxx_std_17)

add_executable(skyspaces_euler3d_reference_demo
    "${SKYSPACES_EULER3D_DEMO_DIR}/Euler3DReferenceDemo.cpp"
)

target_link_libraries(skyspaces_euler3d_reference_demo
    PRIVATE
        skyspaces_euler3d_reference
)

target_compile_features(skyspaces_euler3d_reference_demo PRIVATE cxx_std_17)

if(MSVC)
    target_compile_options(skyspaces_euler3d_reference PRIVATE /W4 /permissive- /utf-8)
    target_compile_options(skyspaces_euler3d_reference_demo PRIVATE /utf-8)
else()
    target_compile_options(skyspaces_euler3d_reference PRIVATE -Wall -Wextra -Wpedantic)
endif()

set_target_properties(skyspaces_euler3d_reference PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
