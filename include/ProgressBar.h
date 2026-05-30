#pragma once

#include <chrono>
#include <string>

namespace skyspaces {

void PrintProgressBar(
    int completed,
    int total,
    const std::string& label,
    const std::string& status = {},
    std::chrono::steady_clock::time_point start_time = {},
    int bar_width = 40);

}  // namespace skyspaces
