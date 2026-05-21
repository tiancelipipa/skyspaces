#pragma once

#include "MathTypes.h"

namespace skyspaces {

enum class InterpolationMethod2D {
    Nearest,
    Bilinear,
    BicubicCatmullRom,
    BicubicLagrange,
};

// Interpolation utilities for 2D scalar grids.
// Provides several sampling methods and short notes on trade-offs.
// Sample a scalar grid stored as an Eigen row-major array (width x height).
// Coordinates (x,y) are in grid index space: 0..width-1, 0..height-1.
Real Sample(
    const ScalarArray2D& data,
    int width,
    int height,
    Real x,
    Real y,
    InterpolationMethod2D method);
Real SampleNearest(const ScalarArray2D& data, int width, int height, Real x, Real y);
Real SampleBilinear(const ScalarArray2D& data, int width, int height, Real x, Real y);
Real SampleBicubicCatmullRom(const ScalarArray2D& data, int width, int height, Real x, Real y);
Real SampleBicubicLagrange(const ScalarArray2D& data, int width, int height, Real x, Real y);

} // namespace skyspaces
