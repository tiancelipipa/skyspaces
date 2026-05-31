@echo off
setlocal EnableDelayedExpansion

set "DEMO=%~1"
set "BUILD_CONFIG=%~2"

if "%DEMO%"=="" set "DEMO=smoke_density"

rem Keep the old usage working: run.bat Debug or run.bat Release.
if /I "%DEMO%"=="Release" (
    set "BUILD_CONFIG=Release"
    set "DEMO=smoke_density"
)
if /I "%DEMO%"=="Debug" (
    set "BUILD_CONFIG=Debug"
    set "DEMO=smoke_density"
)

if "%BUILD_CONFIG%"=="" set "BUILD_CONFIG=Release"

if /I not "%BUILD_CONFIG%"=="Release" if /I not "%BUILD_CONFIG%"=="Debug" (
    call :usage
    exit /b 1
)

set "TARGET="
set "EXE="
set "EXE_PATH="
set "ARGS="

if /I "%DEMO%"=="smoke_density" (
    set "TARGET=skyspaces_smoke_density_demo"
    set "EXE=skyspaces_smoke_density_demo.exe"
    set "ARGS=outputs\smoke_density\smoke_frames 1800 512 512 outputs\smoke_density\smoke_density.mp4 60"
) else if /I "%DEMO%"=="inter_layer_convection" (
    set "TARGET=skyspaces_inter_layer_convection_demo"
    set "EXE=skyspaces_inter_layer_convection_demo.exe"
    set "ARGS=outputs\inter_layer_convection\frames 96"
) else if /I "%DEMO%"=="cloud_layer_volume" (
    set "TARGET=skyspaces_cloud_layer_volume_demo"
    set "EXE=skyspaces_cloud_layer_volume_demo.exe"
    set "ARGS=outputs\cloud_layer_animation 120"
) else if /I "%DEMO%"=="euler3d_reference" (
    set "TARGET=skyspaces_euler3d_reference_demo"
    set "EXE=skyspaces_euler3d_reference_demo.exe"
    if "%EULER3D_RES%"=="" set "EULER3D_RES=128"
    if "%EULER3D_FRAMES%"=="" set "EULER3D_FRAMES=64"
    if "%EULER3D_WARMUP%"=="" set "EULER3D_WARMUP=0"
    set "ARGS=!EULER3D_RES! !EULER3D_RES! !EULER3D_RES! !EULER3D_FRAMES! !EULER3D_WARMUP! outputs\euler3d"
) else (
    echo Unknown demo: %DEMO%
    call :usage
    exit /b 1
)

if not exist build\CMakeCache.txt (
    cmake -S . -B build
    if errorlevel 1 exit /b %errorlevel%
)

cmake --build build --target %TARGET% --config "%BUILD_CONFIG%"
if errorlevel 1 exit /b %errorlevel%

set "EXE_PATH=.\build\%BUILD_CONFIG%\%EXE%"
if /I "%DEMO%"=="euler3d_reference" (
    set "EXE_PATH=.\build\reference\euler3d\%BUILD_CONFIG%\%EXE%"
)

"%EXE_PATH%" %ARGS%
exit /b %errorlevel%

:usage
echo Usage: run.bat [demo] [Release^|Debug]
echo.
echo Demos:
echo   smoke_density            Build and run skyspaces_smoke_density_demo ^(default^)
echo   inter_layer_convection   Build and run skyspaces_inter_layer_convection_demo
echo   cloud_layer_volume       Build and run skyspaces_cloud_layer_volume_demo
echo   euler3d_reference        Build and run skyspaces_euler3d_reference_demo and profiling plots
echo.
echo Examples:
echo   run.bat
echo   run.bat Debug
echo   run.bat smoke_density Debug
echo   run.bat euler3d_reference Release
echo.
echo Euler3D profiling can be sized with environment variables:
echo   set EULER3D_RES=48
echo   set EULER3D_FRAMES=24
echo   set EULER3D_WARMUP=4
exit /b 0
