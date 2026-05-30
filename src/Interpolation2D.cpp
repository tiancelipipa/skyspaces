#include "Interpolation2D.h"

#include <algorithm>
#include <cmath>

namespace skyspaces {

namespace {

// Sample a scalar grid stored as an Eigen array with rows as x and columns as y.
inline Real GetAt(const ScalarArray2D& d, int w, int h, int x, int y) {
    x = std::clamp(x, 0, w - 1);
    y = std::clamp(y, 0, h - 1);
    return d(x, y);
}

// Catmull-Rom cubic interpolation (1D)
inline Real CatmullRom(Real p0, Real p1, Real p2, Real p3, Real t) {
    Real a = -0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3;
    Real b = p0 - 2.5*p1 + 2.0*p2 - 0.5*p3;
    Real c = -0.5*p0 + 0.5*p2;
    Real d = p1;
    return ((a*t + b)*t + c)*t + d;
}

// Lagrange cubic weight for four-point interpolation at parameter s in [0,1].
// Uses points q_{i-1}, q_i, q_{i+1}, q_{i+2}.
inline void LagrangeWeights(Real s, Real weights[4]) {
    Real s2 = s * s;
    Real s3 = s2 * s;
    weights[0] = -s3 / 6.0 + s2 / 2.0 - s / 3.0;
    weights[1] = 1.0 - 0.5 * s - s2 + 0.5 * s3;
    weights[2] = s + 0.5 * s2 - 0.5 * s3;
    weights[3] = s3 / 6.0 - s / 6.0;
}

} // anonymous

Real Sample(
    const ScalarArray2D& data,
    int width,
    int height,
    Real x,
    Real y,
    InterpolationMethod2D method) {
    switch (method) {
        case InterpolationMethod2D::Nearest:
            return SampleNearest(data, width, height, x, y);
        case InterpolationMethod2D::Bilinear:
            return SampleBilinear(data, width, height, x, y);
        case InterpolationMethod2D::BicubicCatmullRom:
            return SampleBicubicCatmullRom(data, width, height, x, y);
        case InterpolationMethod2D::BicubicLagrange:
            return SampleBicubicLagrange(data, width, height, x, y);
    }

    return SampleBilinear(data, width, height, x, y);
}

Real SampleNearest(const ScalarArray2D& data, int width, int height, Real x, Real y) {
    int xi = static_cast<int>(std::round(std::clamp(x, 0.0, static_cast<Real>(width - 1))));
    int yi = static_cast<int>(std::round(std::clamp(y, 0.0, static_cast<Real>(height - 1))));
    return GetAt(data, width, height, xi, yi);
}

Real SampleBilinear(const ScalarArray2D& data, int width, int height, Real x, Real y) {
    // Clamp-to-edge sampling keeps advection stable near solid boundaries.
    Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    int x0 = static_cast<int>(std::floor(cx));
    int y0 = static_cast<int>(std::floor(cy));
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    Real tx = cx - x0;
    Real ty = cy - y0;

    Real c00 = GetAt(data, width, height, x0, y0);
    Real c10 = GetAt(data, width, height, x1, y0);
    Real c01 = GetAt(data, width, height, x0, y1);
    Real c11 = GetAt(data, width, height, x1, y1);

    Real c0 = (1.0 - tx) * c00 + tx * c10;
    Real c1 = (1.0 - tx) * c01 + tx * c11;
    return (1.0 - ty) * c0 + ty * c1;
}

Real SampleBicubicCatmullRom(const ScalarArray2D& data, int width, int height, Real x, Real y) {
    Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    int ix = static_cast<int>(std::floor(cx));
    int iy = static_cast<int>(std::floor(cy));
    Real tx = cx - ix;
    Real ty = cy - iy;

    // Gather a clamped 4x4 neighborhood and interpolate first in x, then in y.
    Real cols[4];
    Real tmp[4];
    for (int j = -1; j <= 2; ++j) {
        int yy = iy + j;
        for (int i = -1; i <= 2; ++i) {
            cols[i+1] = GetAt(data, width, height, ix + i, yy);
        }
        // interp in x for this row
        tmp[j+1] = CatmullRom(cols[0], cols[1], cols[2], cols[3], tx);
    }

    // interp in y using results
    return CatmullRom(tmp[0], tmp[1], tmp[2], tmp[3], ty);
}

Real SampleBicubicLagrange(const ScalarArray2D& data, int width, int height, Real x, Real y) {
    Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    int ix = static_cast<int>(std::floor(cx));
    int iy = static_cast<int>(std::floor(cy));
    Real tx = cx - ix;
    Real ty = cy - iy;

    Real wx[4];
    Real wy[4];
    LagrangeWeights(tx, wx);
    LagrangeWeights(ty, wy);

    Real row[4];
    for (int j = 0; j < 4; ++j) {
        int yy = iy + j - 1;
        Real v0 = GetAt(data, width, height, ix - 1, yy);
        Real v1 = GetAt(data, width, height, ix, yy);
        Real v2 = GetAt(data, width, height, ix + 1, yy);
        Real v3 = GetAt(data, width, height, ix + 2, yy);
        row[j] = wx[0] * v0 + wx[1] * v1 + wx[2] * v2 + wx[3] * v3;
    }

    return wy[0] * row[0] + wy[1] * row[1] + wy[2] * row[2] + wy[3] * row[3];
}

} // namespace skyspaces
