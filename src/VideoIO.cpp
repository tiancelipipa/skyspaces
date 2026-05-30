#include "VideoIO.h"

#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace skyspaces {

namespace {

bool IsFfmpegAvailable() {
#ifdef _WIN32 
/************************************************************
ffmpeg -version
       │
       ├─→ 标准输出(1)  → [屏幕] 
       └─→ 标准错误(2)  → [屏幕]

加上 > NUL 后：
       ├─→ 标准输出(1)  → NUL (丢弃)
       └─→ 标准错误(2)  → [屏幕] ❌ 错误信息还在！

加上 2>&1 后：
       ├─→ 标准错误(2)  ──&1──→ 标准输出(1) 的位置
                                          │
                                          ↓
                                        NUL (一起丢弃)

“> NUL：将标准输出（正常信息）重定向到“空设备”（丢弃）；
2>&1：将标准错误（报错信息）也重定向到标准输出，最终也进入 NUL；
&：表示后面跟的不是文件名，而是文件描述符编号”

************************************************************/
return std::system("ffmpeg -version > NUL 2>&1") == 0;
#else
    return std::system("ffmpeg -version > /dev/null 2>&1") == 0;
#endif
}

}

std::string QuoteCommandArgument(const std::filesystem::path& path) {
    std::string argument = path.generic_string();
    std::string quoted = "\"";
    for (char ch : argument) {
        if (ch == '"') {
            quoted += "\\\"";   // 把路径内的 " 转义成 \"
        } else {
            quoted += ch;
        }
    }
    quoted += "\"";
    return quoted;
}

bool EncodePngSequenceWithFfmpeg(
    const std::filesystem::path& frame_pattern,
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
            << " -y"    // 覆盖旧文件
            << " -loglevel error"
            << " -framerate " << frame_rate
            << " -i " << QuoteCommandArgument(frame_pattern)
            << " -c:v libx264"      // 视频编码器（codec）使用 libx264，用 H.264 标准来压缩视频
            << " -pix_fmt yuv420p"  // 指定视频中的像素格式，yuv420p：色彩空间 YUV，最通用的格式，兼容性极佳
            << " " << QuoteCommandArgument(video_path);

    return std::system(command.str().c_str()) == 0;
}

}  // namespace skyspaces
