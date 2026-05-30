#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace skyspaces {

uint8_t ToByte(double value);

std::filesystem::path NumberedFramePath(
    const std::filesystem::path& output_dir,
    const std::string& prefix,
    int frame,
    const std::string& extension);

std::filesystem::path NumberedFramePattern(
    const std::filesystem::path& output_dir,
    const std::string& prefix,
    const std::string& extension);

void WritePngRgba(
    const std::filesystem::path& path,
    int width,
    int height,
    const std::vector<uint8_t>& rgba);

}  // namespace skyspaces
