#include "Fluid2D.h"
#include "ImageIO.h"
#include "ProgressBar.h"
#include "VideoIO.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kDefaultImageWidth = 512;
constexpr int kDefaultImageHeight = 512;
constexpr int kDefaultFrameCount = 60;
constexpr int kDefaultVideoFrameRate = 60;
constexpr const char* kFramePrefix = "smoke_density";
constexpr const char* kDefaultFramePath = "smoke_frames";
constexpr const char* kDefaultVideoName = "smoke_density.mp4";

struct DemoArguments {
    std::filesystem::path output_dir;
    int frame_count = kDefaultFrameCount;
    int image_width = kDefaultImageWidth;
    int image_height = kDefaultImageHeight;
    std::filesystem::path video_path;
    int video_frame_rate = kDefaultVideoFrameRate;
};

inline int PositiveArgumentOrDefault(int argc, char** argv, int index, int default_value) {
    return argc > index ? std::max(1, std::atoi(argv[index])) : default_value;
}

// 命令行参数：frames_output_dir, frame_count, image_width, image_height, video_path, video_frame_rate
DemoArguments ParseCommandLineArguments(int argc, char** argv) {
    const std::filesystem::path output_dir = argc > 1
        ? std::filesystem::path(argv[1])
        : std::filesystem::current_path() / kDefaultFramePath;

    return {
        output_dir,
        PositiveArgumentOrDefault(argc, argv, 2, kDefaultFrameCount),
        PositiveArgumentOrDefault(argc, argv, 3, kDefaultImageWidth),
        PositiveArgumentOrDefault(argc, argv, 4, kDefaultImageHeight),
        argc > 5 ? std::filesystem::path(argv[5]) : output_dir / kDefaultVideoName,
        PositiveArgumentOrDefault(argc, argv, 6, kDefaultVideoFrameRate)
    };
}

std::array<uint8_t, 4> DensityColor(double density) {
    const double d = std::clamp(density, 0.0, 1.0);
    const double glow = std::clamp((d - 0.15) / 0.85, 0.0, 1.0);

    const uint8_t r = skyspaces::ToByte(0.05 + 0.95 * glow);
    const uint8_t g = skyspaces::ToByte(0.07 + 0.72 * std::sqrt(d));
    const uint8_t b = skyspaces::ToByte(0.10 + 0.45 * d);
    const uint8_t a = 255;
    return {r, g, b, a};
}

skyspaces::CellCenteredScalarGrid2D CreateIrregularCenterSolidLevelSet(
    int resolution_x,
    int resolution_y,
    double cell_size,
    double center_x,
    double center_y) {
    skyspaces::CellCenteredScalarGrid2D phi(resolution_x, resolution_y, 1.0);
    const double base_radius = 0.135;

    for (int x = 0; x < resolution_x; ++x) {
        for (int y = 0; y < resolution_y; ++y) {
            const double px = (static_cast<double>(x) + 0.5) * cell_size;
            const double py = (static_cast<double>(y) + 0.5) * cell_size;
            const double dx = px - center_x;
            const double dy = py - center_y;
            const double angle = std::atan2(dy, dx);
            const double radius =
                base_radius *
                (1.0 +
                 0.24 * std::sin(5.0 * angle + 0.35) +
                 0.13 * std::sin(9.0 * angle - 0.80) +
                 0.07 * std::cos(13.0 * angle + 1.40));
            phi(x, y) = std::sqrt(dx * dx + dy * dy) - radius;
        }
    }

    return phi;
}

skyspaces::CellCenteredScalarGrid2D CreateAnimatedIrregularSolidLevelSet(
    int resolution_x,
    int resolution_y,
    double cell_size,
    double time) {
    const double center_x = 0.52 + 0.070 * std::sin(0.95 * time);
    const double center_y = 0.50 + 0.060 * std::sin(1.35 * time + 0.75);
    return CreateIrregularCenterSolidLevelSet(resolution_x, resolution_y, cell_size, center_x, center_y);
}

skyspaces::FaceCenteredVectorGrid2D CreateAnimatedSolidVelocity(
    int resolution_x,
    int resolution_y,
    double cell_size,
    double time) {
    skyspaces::FaceCenteredVectorGrid2D solid_velocity(resolution_x, resolution_y, 0.0);
    const double center_velocity_x = 0.070 * 0.95 * std::cos(0.95 * time);
    const double center_velocity_y = 0.060 * 1.35 * std::cos(1.35 * time + 0.75);

    for (int x = 0; x < solid_velocity.UWidth(); ++x) {
        for (int y = 0; y < solid_velocity.UHeight(); ++y) {
            solid_velocity.U(x, y) = center_velocity_x;
        }
    }
    for (int x = 0; x < solid_velocity.VWidth(); ++x) {
        for (int y = 0; y < solid_velocity.VHeight(); ++y) {
            solid_velocity.V(x, y) = center_velocity_y;
        }
    }

    return solid_velocity;
}

std::vector<uint8_t> RenderDensityImage(
    const skyspaces::Fluid2D& fluid,
    double max_smoke_density,
    int image_width,
    int image_height) {
    std::vector<uint8_t> pixels(static_cast<std::size_t>(image_width) *
                                static_cast<std::size_t>(image_height) * 4u);
    const skyspaces::CellCenteredScalarGrid2D& solid_phi = fluid.SolidLevelSet();
    const bool has_level_set_solid =
        fluid.Solid().Mode() == skyspaces::SolidBoundaryMode2D::LevelSet;
    const double solid_outline_half_width = 0.35 * fluid.CellSize();

    for (int y = 0; y < image_height; ++y) {
        for (int x = 0; x < image_width; ++x) {
            const double u = (static_cast<double>(x) + 0.5) / static_cast<double>(image_width);
            const double v = 1.0 - (static_cast<double>(y) + 0.5) / static_cast<double>(image_height);
            std::array<uint8_t, 4> color = DensityColor(fluid.SampleSmokeNormalized(u, v) / max_smoke_density);

            if (has_level_set_solid) {
                const int sx = std::clamp(
                    static_cast<int>(u * static_cast<double>(solid_phi.Width())),
                    0,
                    solid_phi.Width() - 1);
                const int sy = std::clamp(
                    static_cast<int>(v * static_cast<double>(solid_phi.Height())),
                    0,
                    solid_phi.Height() - 1);
                const double phi = solid_phi(sx, sy);
                if (std::abs(phi) < solid_outline_half_width) {
                    color = {185, 195, 205, 255};
                }
            }

            const std::size_t offset =
                (static_cast<std::size_t>(y) * static_cast<std::size_t>(image_width) +
                 static_cast<std::size_t>(x)) *
                4u;
            pixels[offset + 0] = color[0];
            pixels[offset + 1] = color[1];
            pixels[offset + 2] = color[2];
            pixels[offset + 3] = color[3];
        }
    }

    return pixels;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        skyspaces::Fluid2D::Config config;
        config.numeric_epsilon = 1e-12;
        config.resolution_x = 48;
        config.resolution_y = 48;
        config.cell_size = 1.0 / 48.0;
        config.time_step = 1.0 / 60.0;
        config.max_substeps = 8;
        config.cfl_number = 1.0;

        config.fluid_density = 1.0;
        config.pressure_iterations = 400;
        config.pressure_tolerance = 5e-6;

        config.ambient_temperature = 273.15;
        config.smoke_dissipation = 0.015;
        config.temperature_dissipation = 0.025;
        config.max_smoke_density = 1.0;

        config.source_min_x = 0.04;
        config.source_max_x = 0.10;
        config.source_min_y = 0.44;
        config.source_max_y = 0.56;
        config.source_smoke_rate = 3.0;
        config.source_temperature = 300.0;
        config.source_temperature_rate = 1.0;
        config.source_acceleration_x = 2.0;
        config.source_acceleration_y = 0.0;

        config.vorticity_confinement = 0.0; // 无粘流体

        config.advection_integrator = skyspaces::AdvectionIntegrator2D::ExplicitEuler;
        config.advection_interpolation = skyspaces::InterpolationMethod2D::Bilinear;
        config.advection_scheme = skyspaces::AdvectionScheme2D::MacCormackBFECC;
        config.pressure_solver = skyspaces::PressureSolver2D::ConjugateGradient;

        skyspaces::Fluid2D fluid(config);

        const DemoArguments args = ParseCommandLineArguments(argc, argv);
        const int frame_count = args.frame_count;
        const int image_width = args.image_width;
        const int image_height = args.image_height;
        const std::filesystem::path& output_dir = args.output_dir;
        const std::filesystem::path& video_path = args.video_path;
        const int video_frame_rate = args.video_frame_rate;
        std::filesystem::create_directories(output_dir);

        skyspaces::PrintProgressBar(0, frame_count, "Writing PNG frames");
        for (int frame = 0; frame < frame_count; ++frame) {
            fluid.SetSolidLevelSet(CreateAnimatedIrregularSolidLevelSet(
                config.resolution_x,
                config.resolution_y,
                config.cell_size,
                fluid.Time()));
            fluid.SetSolidVelocity(CreateAnimatedSolidVelocity(
                config.resolution_x,
                config.resolution_y,
                config.cell_size,
                fluid.Time()));
            fluid.Step();
            skyspaces::WritePngRgba(
                skyspaces::NumberedFramePath(output_dir, kFramePrefix, frame, ".png"),
                image_width,
                image_height,
                RenderDensityImage(
                    fluid,
                    config.max_smoke_density,
                    image_width,
                    image_height));
            skyspaces::PrintProgressBar(frame + 1, frame_count, "Writing PNG frames");
        }
        std::cout << "\n";

        const std::filesystem::path frame_pattern =
            skyspaces::NumberedFramePattern(output_dir, kFramePrefix, ".png");
        const bool video_saved =
            skyspaces::EncodePngSequenceWithFfmpeg(frame_pattern, video_path, video_frame_rate);

        std::cout << "Saved smoke density PNG sequence: " << output_dir.string() << "\n";
        if (video_saved) {
            std::cout << "Saved smoke density video: " << video_path.string() << "\n";
        } else {
            std::cout << "Could not encode video. Install ffmpeg or add it to PATH, then run:\n"
                      << "ffmpeg -y -framerate " << video_frame_rate
                      << " -i " << skyspaces::QuoteCommandArgument(frame_pattern)
                      << " -c:v libx264 -pix_fmt yuv420p "
                      << skyspaces::QuoteCommandArgument(video_path) << "\n";
        }
        std::cout << "Frame count: " << frame_count
                  << ", image size: " << image_width << "x" << image_height
                  << ", time: " << fluid.Time()
                  << ", max smoke density: " << fluid.SmokeDensity().Data().maxCoeff()
                  << ", max velocity: "
                  << std::max(
                         fluid.Velocity().UData().abs().maxCoeff(),
                         fluid.Velocity().VData().abs().maxCoeff())
                  << ", pressure residual: " << fluid.LastPressureResidual() << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Smoke density demo failed: " << ex.what() << "\n";
        return 1;
    }
}
