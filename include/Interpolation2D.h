#pragma once

#include <vector>

namespace skyspaces {

// Interpolation utilities for 2D scalar grids.
// Provides several sampling methods and short notes on trade-offs.
// Sample a scalar grid stored in row-major order (width x height).
// Coordinates (x,y) are in grid index space: 0..width-1, 0..height-1.
double SampleNearest(const std::vector<double>& data, int width, int height, double x, double y);
double SampleBilinear(const std::vector<double>& data, int width, int height, double x, double y);
double SampleBicubicCatmullRom(const std::vector<double>& data, int width, int height, double x, double y);
double SampleBicubicLagrange(const std::vector<double>& data, int width, int height, double x, double y);

} // namespace skyspaces
