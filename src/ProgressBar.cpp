#include "ProgressBar.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace skyspaces {
namespace {

std::string FormatDuration(double seconds) {
    seconds = std::max(0.0, seconds);
    const int total_seconds = static_cast<int>(seconds + 0.5);
    const int minutes = total_seconds / 60;
    const int remaining_seconds = total_seconds % 60;

    std::ostringstream stream;
    stream << minutes << "m " << std::setw(2) << std::setfill('0') << remaining_seconds << "s";
    return stream.str();
}

}  // namespace

void PrintProgressBar(
    int completed,
    int total,
    const std::string& label,
    const std::string& status,
    std::chrono::steady_clock::time_point start_time,
    int bar_width) {
    const int safe_total = std::max(1, total);
    const int safe_bar_width = std::max(1, bar_width);
    const int clamped_completed = std::clamp(completed, 0, safe_total);
    const double fraction =
        static_cast<double>(clamped_completed) / static_cast<double>(safe_total);
    const int filled = static_cast<int>(fraction * static_cast<double>(safe_bar_width) + 0.5);
    const int percent = static_cast<int>(fraction * 100.0 + 0.5);

    std::ostringstream stream;
    stream << "\r";     // 将光标移动到当前行的开头
    if (!label.empty()) {
        stream << label << " ";
    }

    stream << "[";
    for (int i = 0; i < safe_bar_width; ++i) {
        stream << (i < filled ? '#' : '-');
    }
    stream << "] " << std::setw(3) << percent << "% "
           << clamped_completed << "/" << safe_total;

    if (start_time != std::chrono::steady_clock::time_point{}) {
        const auto now = std::chrono::steady_clock::now();
        const double elapsed = std::chrono::duration<double>(now - start_time).count();
        const double eta = clamped_completed > 0
            ? elapsed * static_cast<double>(safe_total - clamped_completed) /
                  static_cast<double>(clamped_completed)
            : 0.0;
        stream << " elapsed " << FormatDuration(elapsed)
               << " eta " << FormatDuration(eta);
    }

    if (!status.empty()) {
        stream << "  " << status << "          ";
    }

    std::cout << stream.str() << std::flush;
}

}  // namespace skyspaces
