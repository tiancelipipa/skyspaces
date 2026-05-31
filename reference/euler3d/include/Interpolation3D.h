#pragma once

#include "MathTypes.h"

#include <vector>

namespace skyspaces::reference {

enum class InterpolationMethod3D {
    Nearest,
    Trilinear,
    TricubicCatmullRom,
    TricubicLagrange,
};

Real Sample(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z,
    InterpolationMethod3D method);
Real SampleNearest(const std::vector<Real>& data, int width, int height, int depth, Real x, Real y, Real z);
Real SampleTrilinear(const std::vector<Real>& data, int width, int height, int depth, Real x, Real y, Real z);
Real SampleTricubicCatmullRom(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z);
Real SampleTricubicLagrange(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z);

}  // namespace skyspaces::reference
