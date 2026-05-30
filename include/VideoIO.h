#pragma once

#include <filesystem>
#include <string>

namespace skyspaces {

// 将文件路径转为命令行参数，用双引号包装（避免路径中的空格出错）
std::string QuoteCommandArgument(const std::filesystem::path& path);

bool EncodePngSequenceWithFfmpeg(
    const std::filesystem::path& frame_pattern,
    const std::filesystem::path& video_path,
    int frame_rate);

}  // namespace skyspaces
