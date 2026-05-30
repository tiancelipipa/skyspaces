@echo off
set "BUILD_CONFIG=%~1"
if "%BUILD_CONFIG%"=="" set "BUILD_CONFIG=Release"

if /I not "%BUILD_CONFIG%"=="Release" if /I not "%BUILD_CONFIG%"=="Debug" (
    echo Usage: run.bat [Release^|Debug]
    exit /b 1
)

if not exist build (
    cmake -S . -B build
    if errorlevel 1 exit /b %errorlevel%
)

cmake --build build --target skyspaces_smoke_density_demo --config "%BUILD_CONFIG%"
if errorlevel 1 exit /b %errorlevel%

@REM Command line: output_dir, frame_count, image_width, image_height, video_path, video_frame_rate
".\build\%BUILD_CONFIG%\skyspaces_smoke_density_demo.exe" outputs\smoke_density\smoke_frames 1800 512 512 outputs\smoke_density\smoke_density.mp4 60

@REM cmake --build build --target skyspaces_inter_layer_convection_demo --config Debug
@REM build\Debug\skyspaces_inter_layer_convection_demo.exe build\inter_layer_convection_test 60

@REM cmake --build build --target skyspaces_cloud_layer_volume_demo --config Release
@REM build\Debug\skyspaces_cloud_layer_volume_demo.exe build\cloud_layer_volume
@REM build\Release\skyspaces_cloud_layer_volume_demo.exe outputs\cloud_layer_animation 120
