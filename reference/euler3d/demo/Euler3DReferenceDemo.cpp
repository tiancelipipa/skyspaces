#include "EulerFluid3D.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Timings = skyspaces::reference::EulerFluid3D::StepTimings;

constexpr int kDefaultResolution = 128;
constexpr int kDefaultMeasuredFrames = 64;
constexpr int kDefaultWarmupFrames = 0;
constexpr const char* kDefaultOutputDir = "outputs/euler3d";

int PositiveArgumentOrDefault(int argc, char** argv, int index, int default_value) {
    return argc > index ? std::max(1, std::atoi(argv[index])) : default_value;
}

skyspaces::Real SecondsSince(Clock::time_point start) {
    return std::chrono::duration<skyspaces::Real>(Clock::now() - start).count();
}

skyspaces::Real MaxAbs(const std::vector<skyspaces::Real>& values) {
    skyspaces::Real result = 0.0;
    for (const skyspaces::Real value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

skyspaces::Real MaxValue(const std::vector<skyspaces::Real>& values) {
    return values.empty() ? 0.0 : *std::max_element(values.begin(), values.end());
}

void AddTimings(Timings& total, const Timings& step) {
    total.substep_count += step.substep_count;
    total.total += step.total;
    total.advect_velocity += step.advect_velocity;
    total.velocity_source += step.velocity_source;
    total.boundary_conditions += step.boundary_conditions;
    total.pre_projection_boundary += step.pre_projection_boundary;
    total.projection_boundary += step.projection_boundary;
    total.projection += step.projection;
    total.projection_divergence += step.projection_divergence;
    total.pressure_system_build += step.pressure_system_build;
    total.pressure_solve += step.pressure_solve;
    total.pressure_apply_gradient += step.pressure_apply_gradient;
    total.projection_validation += step.projection_validation;
    total.advect_scalars += step.advect_scalars;
    total.scalar_advection += step.scalar_advection;
    total.scalar_postprocess += step.scalar_postprocess;
    total.scalar_source += step.scalar_source;
}

skyspaces::Real MillisecondsPerFrame(skyspaces::Real seconds, int measured_frames) {
    return 1000.0 * seconds / static_cast<skyspaces::Real>(measured_frames);
}

void PrintTimingRow(
    const std::string& name,
    skyspaces::Real seconds,
    skyspaces::Real baseline_seconds,
    int measured_frames) {
    const skyspaces::Real percent = baseline_seconds > 0.0 ? 100.0 * seconds / baseline_seconds : 0.0;
    std::cout << std::left << std::setw(34) << name
              << std::right << std::setw(12) << std::fixed << std::setprecision(3)
              << MillisecondsPerFrame(seconds, measured_frames)
              << std::setw(11) << std::fixed << std::setprecision(2) << percent << "%\n";
}

void PrintProgress(const char* label, int completed, int total) {
    if (total <= 0) {
        return;
    }

    const double percent = 100.0 * static_cast<double>(completed) / static_cast<double>(total);
    std::cout << "\r" << label << ": " << completed << "/" << total
              << " (" << std::fixed << std::setprecision(1) << percent << "%)" << std::flush;
    if (completed == total) {
        std::cout << "\n";
    }
}

void WriteFrameCsv(const std::filesystem::path& path, const std::vector<Timings>& frames) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not write frame CSV: " + path.string());
    }

    out << "frame,substeps,total_ms,advect_velocity_ms,velocity_source_ms,boundary_conditions_ms,"
           "pre_projection_boundary_ms,projection_boundary_ms,projection_ms,projection_divergence_ms,"
           "pressure_system_build_ms,pressure_solve_ms,pressure_apply_gradient_ms,projection_validation_ms,"
           "advect_scalars_ms,scalar_advection_ms,scalar_postprocess_ms,scalar_source_ms\n";
    for (std::size_t frame = 0; frame < frames.size(); ++frame) {
        const Timings& t = frames[frame];
        out << frame << ','
            << t.substep_count << ','
            << 1000.0 * t.total << ','
            << 1000.0 * t.advect_velocity << ','
            << 1000.0 * t.velocity_source << ','
            << 1000.0 * t.boundary_conditions << ','
            << 1000.0 * t.pre_projection_boundary << ','
            << 1000.0 * t.projection_boundary << ','
            << 1000.0 * t.projection << ','
            << 1000.0 * t.projection_divergence << ','
            << 1000.0 * t.pressure_system_build << ','
            << 1000.0 * t.pressure_solve << ','
            << 1000.0 * t.pressure_apply_gradient << ','
            << 1000.0 * t.projection_validation << ','
            << 1000.0 * t.advect_scalars << ','
            << 1000.0 * t.scalar_advection << ','
            << 1000.0 * t.scalar_postprocess << ','
            << 1000.0 * t.scalar_source << '\n';
    }
}

void WriteSummaryCsv(const std::filesystem::path& path, const Timings& accumulated, int measured_frames) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Could not write summary CSV: " + path.string());
    }

    const auto row = [&](const char* name, skyspaces::Real seconds) {
        const skyspaces::Real ms_per_frame = MillisecondsPerFrame(seconds, measured_frames);
        const skyspaces::Real percent = accumulated.total > 0.0 ? 100.0 * seconds / accumulated.total : 0.0;
        out << name << ',' << 1000.0 * seconds << ',' << ms_per_frame << ',' << percent << '\n';
    };

    out << "stage,total_ms,ms_per_frame,percent_total\n";
    row("total step", accumulated.total);
    row("advect velocity", accumulated.advect_velocity);
    row("velocity source", accumulated.velocity_source);
    row("boundary conditions total", accumulated.boundary_conditions);
    row("pre-projection boundary", accumulated.pre_projection_boundary);
    row("projection total", accumulated.projection);
    row("divergence evaluations", accumulated.projection_divergence);
    row("pressure system build", accumulated.pressure_system_build);
    row("pressure solve", accumulated.pressure_solve);
    row("apply pressure gradient", accumulated.pressure_apply_gradient);
    row("projection boundary", accumulated.projection_boundary);
    row("projection validation", accumulated.projection_validation);
    row("advect scalars total", accumulated.advect_scalars);
    row("scalar advection scheme", accumulated.scalar_advection);
    row("scalar postprocess", accumulated.scalar_postprocess);
    row("scalar source", accumulated.scalar_source);
}

std::string QuoteCommandArgument(const std::filesystem::path& path) {
    std::string text = path.string();
    std::string quoted = "\"";
    for (const char ch : text) {
        if (ch == '"') {
            quoted += "\\\"";
        } else {
            quoted += ch;
        }
    }
    quoted += "\"";
    return quoted;
}

void RunPlotProfile(const std::filesystem::path& frame_csv_path, const std::filesystem::path& output_dir) {
    const std::filesystem::path script_path =
        std::filesystem::current_path() / "reference" / "euler3d" / "scripts" / "plot_profile.py";
    const std::string command =
        "conda run -n skyspaces python " +
        QuoteCommandArgument(script_path) +
        " --input " +
        QuoteCommandArgument(frame_csv_path) +
        " --out " +
        QuoteCommandArgument(output_dir);
    std::cout << "Generating profiling plots...\n" << std::flush;
    const int exit_code = std::system(command.c_str());
    if (exit_code != 0) {
        throw std::runtime_error("plot_profile.py failed with exit code " + std::to_string(exit_code));
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        skyspaces::reference::EulerFluid3D::Config config;
        config.resolution_x = PositiveArgumentOrDefault(argc, argv, 1, kDefaultResolution);
        config.resolution_y = PositiveArgumentOrDefault(argc, argv, 2, kDefaultResolution);
        config.resolution_z = PositiveArgumentOrDefault(argc, argv, 3, kDefaultResolution);
        config.cell_size = 1.0 / static_cast<skyspaces::Real>(config.resolution_x);
        config.time_step = 1.0 / 60.0;
        config.max_substeps = 8;
        config.cfl_number = 1.0;
        config.pressure_iterations = 300;
        config.pressure_tolerance = 1e-5;
        config.smoke_dissipation = 0.015;
        config.temperature_dissipation = 0.025;
        config.source_min_x = 0.44;
        config.source_max_x = 0.56;
        config.source_min_y = 0.05;
        config.source_max_y = 0.12;
        config.source_min_z = 0.44;
        config.source_max_z = 0.56;
        config.source_acceleration_x = 0.0;
        config.source_acceleration_y = 2.0;
        config.source_acceleration_z = 0.0;
        config.advection_integrator = skyspaces::reference::AdvectionIntegrator3D::ExplicitEuler;
        config.advection_interpolation = skyspaces::reference::InterpolationMethod3D::Trilinear;
        config.advection_scheme = skyspaces::reference::AdvectionScheme3D::MacCormackBFECC;
        config.pressure_solver = skyspaces::reference::PressureSolver3D::ConjugateGradient;

        const int measured_frames = PositiveArgumentOrDefault(argc, argv, 4, kDefaultMeasuredFrames);
        const int warmup_frames = PositiveArgumentOrDefault(argc, argv, 5, kDefaultWarmupFrames);
        const std::filesystem::path output_dir =
            argc > 6 ? std::filesystem::path(argv[6]) : std::filesystem::current_path() / kDefaultOutputDir;
        std::filesystem::create_directories(output_dir);
        skyspaces::reference::EulerFluid3D fluid(config);

        std::cout << "3D Euler reference profiling\n"
                  << "resolution: " << config.resolution_x << "x" << config.resolution_y << "x"
                  << config.resolution_z << "\n"
                  << "warmup frames: " << warmup_frames << "\n"
                  << "measured frames: " << measured_frames << "\n"
                  << "output dir: " << output_dir.string() << "\n"
                  << "Starting simulation...\n"
                  << std::flush;

        for (int frame = 0; frame < warmup_frames; ++frame) {
            fluid.Step();
            PrintProgress("Warmup", frame + 1, warmup_frames);
        }

        Timings accumulated;
        std::vector<Timings> measured_frame_timings;
        measured_frame_timings.reserve(static_cast<std::size_t>(measured_frames));
        const auto measured_start = Clock::now();
        for (int frame = 0; frame < measured_frames; ++frame) {
            fluid.Step();
            const Timings& step_timings = fluid.LastStepTimings();
            measured_frame_timings.push_back(step_timings);
            AddTimings(accumulated, step_timings);
            PrintProgress("Measured frames", frame + 1, measured_frames);
        }
        const skyspaces::Real measured_wall_seconds = SecondsSince(measured_start);

        const std::filesystem::path frame_csv_path = output_dir / "euler3d_profile_frames.csv";
        const std::filesystem::path summary_csv_path = output_dir / "euler3d_profile_summary.csv";
        std::cout << "Writing profiling CSV files...\n" << std::flush;
        WriteFrameCsv(frame_csv_path, measured_frame_timings);
        WriteSummaryCsv(summary_csv_path, accumulated, measured_frames);
        RunPlotProfile(frame_csv_path, output_dir);

        const auto& velocity = fluid.Velocity();
        const skyspaces::Real max_velocity =
            std::max({MaxAbs(velocity.UData()), MaxAbs(velocity.VData()), MaxAbs(velocity.WData())});
        const skyspaces::Real avg_substeps =
            static_cast<skyspaces::Real>(accumulated.substep_count) / static_cast<skyspaces::Real>(measured_frames);

        std::cout << "\nProfiling summary\n"
                  << "resolution: " << config.resolution_x << "x" << config.resolution_y << "x"
                  << config.resolution_z << "\n"
                  << "warmup frames: " << warmup_frames << "\n"
                  << "measured frames: " << measured_frames << "\n"
                  << "avg substeps/frame: " << std::fixed << std::setprecision(2) << avg_substeps << "\n"
                  << "simulated time: " << std::fixed << std::setprecision(4) << fluid.Time() << "\n"
                  << "measured wall time: " << std::fixed << std::setprecision(3) << measured_wall_seconds << " s\n"
                  << "instrumented avg frame: " << std::fixed << std::setprecision(3)
                  << MillisecondsPerFrame(accumulated.total, measured_frames) << " ms\n"
                  << "wall avg frame: " << std::fixed << std::setprecision(3)
                  << MillisecondsPerFrame(measured_wall_seconds, measured_frames) << " ms\n"
                  << "max smoke density: " << MaxValue(fluid.SmokeDensity().Data()) << "\n"
                  << "max velocity: " << max_velocity << "\n"
                  << "last pressure iterations: " << fluid.LastPressureIterations() << "\n"
                  << "last pressure residual: " << std::scientific << std::setprecision(6)
                  << fluid.LastPressureResidual() << std::defaultfloat << "\n"
                  << "frame CSV: " << frame_csv_path.string() << "\n"
                  << "summary CSV: " << summary_csv_path.string() << "\n\n";

        std::cout << std::left << std::setw(34) << "stage"
                  << std::right << std::setw(12) << "ms/frame"
                  << std::setw(12) << "% total\n";
        std::cout << std::string(58, '-') << "\n";
        PrintTimingRow("total step", accumulated.total, accumulated.total, measured_frames);
        PrintTimingRow("advect velocity", accumulated.advect_velocity, accumulated.total, measured_frames);
        PrintTimingRow("velocity source", accumulated.velocity_source, accumulated.total, measured_frames);
        PrintTimingRow("boundary conditions total", accumulated.boundary_conditions, accumulated.total, measured_frames);
        PrintTimingRow("  pre-projection boundary", accumulated.pre_projection_boundary, accumulated.total, measured_frames);
        PrintTimingRow("projection total", accumulated.projection, accumulated.total, measured_frames);
        PrintTimingRow("  divergence evaluations", accumulated.projection_divergence, accumulated.total, measured_frames);
        PrintTimingRow("  pressure system build", accumulated.pressure_system_build, accumulated.total, measured_frames);
        PrintTimingRow("  pressure solve", accumulated.pressure_solve, accumulated.total, measured_frames);
        PrintTimingRow("  apply pressure gradient", accumulated.pressure_apply_gradient, accumulated.total, measured_frames);
        PrintTimingRow("  projection boundary", accumulated.projection_boundary, accumulated.total, measured_frames);
        PrintTimingRow("  projection validation", accumulated.projection_validation, accumulated.total, measured_frames);
        PrintTimingRow("advect scalars total", accumulated.advect_scalars, accumulated.total, measured_frames);
        PrintTimingRow("  scalar advection scheme", accumulated.scalar_advection, accumulated.total, measured_frames);
        PrintTimingRow("  scalar postprocess", accumulated.scalar_postprocess, accumulated.total, measured_frames);
        PrintTimingRow("scalar source", accumulated.scalar_source, accumulated.total, measured_frames);

        return 0;
    } catch (const std::bad_alloc&) {
        std::cerr << "3D Euler reference profiling failed: out of memory while allocating 3D grids or pressure solver data.\n"
                  << "Try a smaller resolution, for example:\n"
                  << "  skyspaces_euler3d_reference_demo.exe 32 32 32 24 4 outputs\\euler3d\n";
        return 1;
    } catch (const std::exception& ex) {
        std::cerr << "3D Euler reference profiling failed: " << ex.what() << "\n";
        return 1;
    }
}
