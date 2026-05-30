#include "ImageIO.h"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace skyspaces {

uint8_t ToByte(double value) {
    return static_cast<uint8_t>(std::clamp(value, 0.0, 1.0) * 255.0 + 0.5);
}

std::filesystem::path NumberedFramePath(
    const std::filesystem::path& output_dir,
    const std::string& prefix,
    int frame,
    const std::string& extension) {
    std::ostringstream name;
    // 设置下一个输出字段的最小宽度为 4 个字符，不足部分用 '0' 填充
    name << prefix << "_" << std::setw(4) << std::setfill('0') << frame << extension;
    return output_dir / name.str();
}

std::filesystem::path NumberedFramePattern(
    const std::filesystem::path& output_dir,
    const std::string& prefix,
    const std::string& extension) {
    return output_dir / (prefix + "_%04d" + extension);
}

void WritePngRgba(
    const std::filesystem::path& path,
    int width,
    int height,
    const std::vector<uint8_t>& rgba) {
    if (width <= 0 || height <= 0) {
        throw std::invalid_argument("PNG dimensions must be positive.");
    }
    if (rgba.size() != static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u) {
        throw std::invalid_argument("RGBA buffer size does not match PNG dimensions.");
    }

    const std::string output_path = path.string();
    const int stride_bytes = width * 4;
    const int ok = stbi_write_png(
        output_path.c_str(),    // 输出文件名
        width,                  // 图像宽度（像素）
        height,                 // 图像高度（像素）
        4,                      // 通道数（1=灰度, 2=灰度+Alpha, 3=RGB, 4=RGBA）
        rgba.data(),            // 像素数据指针
        stride_bytes            // 每行字节数（0 = 自动计算）
    );
    if (ok == 0) {
        throw std::runtime_error("Failed to write PNG output path: " + output_path);
    }
}

}  // namespace skyspaces
