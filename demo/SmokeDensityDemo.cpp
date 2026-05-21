#include "Fluid2D.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr int kImageWidth = 512;
constexpr int kImageHeight = 512;

uint32_t Crc32(const uint8_t* data, std::size_t size) {
    uint32_t crc = 0xffffffffu;
    for (std::size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit) {
            crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
        }
    }
    return crc ^ 0xffffffffu;
}

uint32_t Adler32(const std::vector<uint8_t>& data) {
    constexpr uint32_t mod = 65521u;
    uint32_t a = 1u;
    uint32_t b = 0u;
    for (uint8_t value : data) {
        a = (a + value) % mod;
        b = (b + a) % mod;
    }
    return (b << 16) | a;
}

void AppendBigEndian32(std::vector<uint8_t>& out, uint32_t value) {
    out.push_back(static_cast<uint8_t>((value >> 24) & 0xffu));
    out.push_back(static_cast<uint8_t>((value >> 16) & 0xffu));
    out.push_back(static_cast<uint8_t>((value >> 8) & 0xffu));
    out.push_back(static_cast<uint8_t>(value & 0xffu));
}

void AppendLittleEndian16(std::vector<uint8_t>& out, uint16_t value) {
    out.push_back(static_cast<uint8_t>(value & 0xffu));
    out.push_back(static_cast<uint8_t>((value >> 8) & 0xffu));
}

void AppendChunk(std::vector<uint8_t>& png, const std::array<char, 4>& type, const std::vector<uint8_t>& data) {
    AppendBigEndian32(png, static_cast<uint32_t>(data.size()));

    const std::size_t type_begin = png.size();
    for (char ch : type) {
        png.push_back(static_cast<uint8_t>(ch));
    }

    png.insert(png.end(), data.begin(), data.end());

    const uint32_t crc = Crc32(png.data() + type_begin, 4 + data.size());
    AppendBigEndian32(png, crc);
}

std::vector<uint8_t> DeflateNoCompression(const std::vector<uint8_t>& input) {
    std::vector<uint8_t> out;
    out.reserve(input.size() + input.size() / 65535 + 16);

    out.push_back(0x78);
    out.push_back(0x01);

    std::size_t offset = 0;
    while (offset < input.size()) {
        const std::size_t block_size = std::min<std::size_t>(65535, input.size() - offset);
        const bool final_block = offset + block_size == input.size();

        out.push_back(final_block ? 0x01 : 0x00);
        AppendLittleEndian16(out, static_cast<uint16_t>(block_size));
        AppendLittleEndian16(out, static_cast<uint16_t>(~static_cast<uint16_t>(block_size)));
        out.insert(out.end(), input.begin() + static_cast<std::ptrdiff_t>(offset),
                   input.begin() + static_cast<std::ptrdiff_t>(offset + block_size));

        offset += block_size;
    }

    AppendBigEndian32(out, Adler32(input));
    return out;
}

void WritePngRgba(const std::filesystem::path& path, int width, int height, const std::vector<uint8_t>& rgba) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("PNG dimensions must be positive.");
    }
    if (rgba.size() != static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u) {
        throw std::invalid_argument("RGBA buffer size does not match PNG dimensions.");
    }

    std::vector<uint8_t> scanlines;
    scanlines.reserve(static_cast<std::size_t>(height) * (1u + static_cast<std::size_t>(width) * 4u));

    for (int y = 0; y < height; ++y) {
        scanlines.push_back(0);
        const std::size_t row_begin = static_cast<std::size_t>(y) * static_cast<std::size_t>(width) * 4u;
        scanlines.insert(scanlines.end(), rgba.begin() + static_cast<std::ptrdiff_t>(row_begin),
                         rgba.begin() + static_cast<std::ptrdiff_t>(row_begin + static_cast<std::size_t>(width) * 4u));
    }

    std::vector<uint8_t> png;
    png.insert(png.end(), {0x89, 'P', 'N', 'G', '\r', '\n', 0x1a, '\n'});

    std::vector<uint8_t> ihdr;
    AppendBigEndian32(ihdr, static_cast<uint32_t>(width));
    AppendBigEndian32(ihdr, static_cast<uint32_t>(height));
    ihdr.push_back(8);
    ihdr.push_back(6);
    ihdr.push_back(0);
    ihdr.push_back(0);
    ihdr.push_back(0);

    AppendChunk(png, {'I', 'H', 'D', 'R'}, ihdr);
    AppendChunk(png, {'I', 'D', 'A', 'T'}, DeflateNoCompression(scanlines));
    AppendChunk(png, {'I', 'E', 'N', 'D'}, {});

    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open PNG output path: " + path.string());
    }
    file.write(reinterpret_cast<const char*>(png.data()), static_cast<std::streamsize>(png.size()));
}

uint8_t ToByte(double value) {
    const double clamped = std::clamp(value, 0.0, 1.0);
    return static_cast<uint8_t>(clamped * 255.0 + 0.5);
}

std::array<uint8_t, 4> DensityColor(double density) {
    const double d = std::clamp(density, 0.0, 1.0);
    const double glow = std::clamp((d - 0.15) / 0.85, 0.0, 1.0);

    const uint8_t r = ToByte(0.05 + 0.95 * glow);
    const uint8_t g = ToByte(0.07 + 0.72 * std::sqrt(d));
    const uint8_t b = ToByte(0.10 + 0.45 * d);
    const uint8_t a = 255;
    return {r, g, b, a};
}

std::filesystem::path OutputDirectory(int argc, char** argv) {
    if (argc > 1) {
        return argv[1];
    }
    return std::filesystem::current_path() / "smoke_frames";
}

int FrameCount(int argc, char** argv) {
    if (argc > 2) {
        return std::max(1, std::atoi(argv[2]));
    }
    return 60;
}

int VideoFrameRate(int argc, char** argv) {
    if (argc > 4) {
        return std::max(1, std::atoi(argv[4]));
    }
    return 60;
}

std::filesystem::path FramePath(const std::filesystem::path& output_dir, int frame) {
    std::ostringstream name;
    name << "smoke_density_" << std::setw(4) << std::setfill('0') << frame << ".png";
    return output_dir / name.str();
}

std::filesystem::path FramePattern(const std::filesystem::path& output_dir) {
    return output_dir / "smoke_density_%04d.png";
}

std::filesystem::path VideoPath(int argc, char** argv, const std::filesystem::path& output_dir) {
    if (argc > 3) {
        return argv[3];
    }
    return output_dir / "smoke_density.mp4";
}

void PrintProgressBar(int completed, int total) {
    constexpr int bar_width = 40;
    const int safe_total = std::max(1, total);
    const double progress =
        static_cast<double>(std::clamp(completed, 0, safe_total)) / static_cast<double>(safe_total);
    const int filled = static_cast<int>(progress * static_cast<double>(bar_width) + 0.5);
    const int percent = static_cast<int>(progress * 100.0 + 0.5);

    std::cout << "\rWriting PNG frames [";
    for (int i = 0; i < bar_width; ++i) {
        std::cout << (i < filled ? '#' : '-');
    }
    std::cout << "] " << std::setw(3) << percent << "% "
              << std::clamp(completed, 0, safe_total) << "/" << safe_total << std::flush;
}

std::string QuoteCommandArgument(const std::filesystem::path& path) {
    std::string argument = path.generic_string();
    std::string quoted = "\"";
    for (char ch : argument) {
        if (ch == '"') {
            quoted += "\\\"";
        } else {
            quoted += ch;
        }
    }
    quoted += "\"";
    return quoted;
}

bool IsFfmpegAvailable() {
#ifdef _WIN32
    return std::system("ffmpeg -version > NUL 2>&1") == 0;
#else
    return std::system("ffmpeg -version > /dev/null 2>&1") == 0;
#endif
}

bool EncodeVideoWithFfmpeg(
    const std::filesystem::path& output_dir,
    const std::filesystem::path& video_path,
    int frame_rate) {
    if (!IsFfmpegAvailable()) {
        return false;
    }

    const std::filesystem::path parent = video_path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }

    std::ostringstream command;
    command << "ffmpeg"
            << " -y"
            << " -loglevel error"
            << " -framerate " << frame_rate
            << " -i " << QuoteCommandArgument(FramePattern(output_dir))
            << " -c:v libx264"
            << " -pix_fmt yuv420p"
            << " " << QuoteCommandArgument(video_path);

    return std::system(command.str().c_str()) == 0;
}

std::vector<uint8_t> RenderDensityImage(const skyspaces::Fluid2D& fluid, double max_smoke_density) {
    std::vector<uint8_t> pixels(static_cast<std::size_t>(kImageWidth) *
                                static_cast<std::size_t>(kImageHeight) * 4u);

    for (int y = 0; y < kImageHeight; ++y) {
        for (int x = 0; x < kImageWidth; ++x) {
            const double u = (static_cast<double>(x) + 0.5) / static_cast<double>(kImageWidth);
            const double v = 1.0 - (static_cast<double>(y) + 0.5) / static_cast<double>(kImageHeight);
            const double density = fluid.SampleSmokeNormalized(u, v) / max_smoke_density;
            const std::array<uint8_t, 4> color = DensityColor(density);

            const std::size_t offset =
                (static_cast<std::size_t>(y) * static_cast<std::size_t>(kImageWidth) +
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
        config.resolution_x = 48;
        config.resolution_y = 48;
        config.cell_size = 1.0 / 48.0;
        config.time_step = 1.0 / 60.0;
        // config.pressure_solver = skyspaces::PressureSolver2D::ConjugateGradient;
        config.pressure_iterations = 400;
        config.pressure_tolerance = 5e-6;
        config.vorticity_confinement = 0.0;
        config.smoke_dissipation = 0.015;
        config.temperature_dissipation = 0.025;
        config.source_min_x = 0.44;
        config.source_max_x = 0.56;
        config.source_min_y = 0.04;
        config.source_max_y = 0.10;
        config.source_smoke_rate = 3.0;
        config.source_temperature = 700.0;
        config.source_temperature_rate = 8.0;
        config.source_acceleration_y = 2.0;

        skyspaces::Fluid2D fluid(config);
        const int frame_count = FrameCount(argc, argv);
        const std::filesystem::path output_dir = OutputDirectory(argc, argv);
        const std::filesystem::path video_path = VideoPath(argc, argv, output_dir);
        const int video_frame_rate = VideoFrameRate(argc, argv);
        std::filesystem::create_directories(output_dir);

        PrintProgressBar(0, frame_count);
        for (int frame = 0; frame < frame_count; ++frame) {
            fluid.Step();
            WritePngRgba(FramePath(output_dir, frame), kImageWidth, kImageHeight,
                         RenderDensityImage(fluid, config.max_smoke_density));
            PrintProgressBar(frame + 1, frame_count);
        }
        std::cout << "\n";

        const bool video_saved = EncodeVideoWithFfmpeg(output_dir, video_path, video_frame_rate);

        std::cout << "Saved smoke density PNG sequence: " << output_dir.string() << "\n";
        if (video_saved) {
            std::cout << "Saved smoke density video: " << video_path.string() << "\n";
        } else {
            std::cout << "Could not encode video. Install ffmpeg or add it to PATH, then run:\n"
                      << "ffmpeg -y -framerate " << video_frame_rate
                      << " -i " << QuoteCommandArgument(FramePattern(output_dir))
                      << " -c:v libx264 -pix_fmt yuv420p "
                      << QuoteCommandArgument(video_path) << "\n";
        }
        std::cout << "Frame count: " << frame_count
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
