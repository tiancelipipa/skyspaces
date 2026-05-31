#include "Interpolation3D.h"

#include <algorithm>
#include <cmath>

namespace skyspaces::reference {
namespace {

inline int ClampIndex(int index, int size) {
    return std::clamp(index, 0, size - 1);
}

inline std::size_t Index3D(int x, int y, int z, int height, int depth) {
    return (static_cast<std::size_t>(x) * static_cast<std::size_t>(height) +
            static_cast<std::size_t>(y)) *
               static_cast<std::size_t>(depth) +
           static_cast<std::size_t>(z);
}

inline Real GetAt(const std::vector<Real>& data, int width, int height, int depth, int x, int y, int z) {
    x = ClampIndex(x, width);
    y = ClampIndex(y, height);
    z = ClampIndex(z, depth);
    return data[Index3D(x, y, z, height, depth)];
}

inline Real CatmullRom(Real p0, Real p1, Real p2, Real p3, Real t) {
    const Real a = -0.5 * p0 + 1.5 * p1 - 1.5 * p2 + 0.5 * p3;
    const Real b = p0 - 2.5 * p1 + 2.0 * p2 - 0.5 * p3;
    const Real c = -0.5 * p0 + 0.5 * p2;
    return ((a * t + b) * t + c) * t + p1;
}

inline void LagrangeWeights(Real s, Real weights[4]) {
    const Real s2 = s * s;
    const Real s3 = s2 * s;
    weights[0] = -s3 / 6.0 + s2 / 2.0 - s / 3.0;
    weights[1] = 1.0 - 0.5 * s - s2 + 0.5 * s3;
    weights[2] = s + 0.5 * s2 - 0.5 * s3;
    weights[3] = s3 / 6.0 - s / 6.0;
}

}  // namespace

Real Sample(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z,
    InterpolationMethod3D method) {
    switch (method) {
        case InterpolationMethod3D::Nearest:
            return SampleNearest(data, width, height, depth, x, y, z);
        case InterpolationMethod3D::Trilinear:
            return SampleTrilinear(data, width, height, depth, x, y, z);
        case InterpolationMethod3D::TricubicCatmullRom:
            return SampleTricubicCatmullRom(data, width, height, depth, x, y, z);
        case InterpolationMethod3D::TricubicLagrange:
            return SampleTricubicLagrange(data, width, height, depth, x, y, z);
    }

    return SampleTrilinear(data, width, height, depth, x, y, z);
}

Real SampleNearest(const std::vector<Real>& data, int width, int height, int depth, Real x, Real y, Real z) {
    const int xi = static_cast<int>(std::round(std::clamp(x, 0.0, static_cast<Real>(width - 1))));
    const int yi = static_cast<int>(std::round(std::clamp(y, 0.0, static_cast<Real>(height - 1))));
    const int zi = static_cast<int>(std::round(std::clamp(z, 0.0, static_cast<Real>(depth - 1))));
    return GetAt(data, width, height, depth, xi, yi, zi);
}

Real SampleTrilinear(const std::vector<Real>& data, int width, int height, int depth, Real x, Real y, Real z) {
    const Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    const Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    const Real cz = std::clamp(z, 0.0, static_cast<Real>(depth - 1));
    const int x0 = static_cast<int>(std::floor(cx));
    const int y0 = static_cast<int>(std::floor(cy));
    const int z0 = static_cast<int>(std::floor(cz));
    const int x1 = ClampIndex(x0 + 1, width);
    const int y1 = ClampIndex(y0 + 1, height);
    const int z1 = ClampIndex(z0 + 1, depth);
    const Real tx = cx - static_cast<Real>(x0);
    const Real ty = cy - static_cast<Real>(y0);
    const Real tz = cz - static_cast<Real>(z0);

    const Real c00 = (1.0 - tx) * GetAt(data, width, height, depth, x0, y0, z0) +
                     tx * GetAt(data, width, height, depth, x1, y0, z0);
    const Real c10 = (1.0 - tx) * GetAt(data, width, height, depth, x0, y1, z0) +
                     tx * GetAt(data, width, height, depth, x1, y1, z0);
    const Real c01 = (1.0 - tx) * GetAt(data, width, height, depth, x0, y0, z1) +
                     tx * GetAt(data, width, height, depth, x1, y0, z1);
    const Real c11 = (1.0 - tx) * GetAt(data, width, height, depth, x0, y1, z1) +
                     tx * GetAt(data, width, height, depth, x1, y1, z1);
    const Real c0 = (1.0 - ty) * c00 + ty * c10;
    const Real c1 = (1.0 - ty) * c01 + ty * c11;
    return (1.0 - tz) * c0 + tz * c1;
}

Real SampleTricubicCatmullRom(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z) {
    const Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    const Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    const Real cz = std::clamp(z, 0.0, static_cast<Real>(depth - 1));
    const int ix = static_cast<int>(std::floor(cx));
    const int iy = static_cast<int>(std::floor(cy));
    const int iz = static_cast<int>(std::floor(cz));
    const Real tx = cx - static_cast<Real>(ix);
    const Real ty = cy - static_cast<Real>(iy);
    const Real tz = cz - static_cast<Real>(iz);

    Real z_rows[4];
    for (int kz = -1; kz <= 2; ++kz) {
        Real y_rows[4];
        for (int jy = -1; jy <= 2; ++jy) {
            Real x_values[4];
            for (int ix_offset = -1; ix_offset <= 2; ++ix_offset) {
                x_values[ix_offset + 1] =
                    GetAt(data, width, height, depth, ix + ix_offset, iy + jy, iz + kz);
            }
            y_rows[jy + 1] = CatmullRom(x_values[0], x_values[1], x_values[2], x_values[3], tx);
        }
        z_rows[kz + 1] = CatmullRom(y_rows[0], y_rows[1], y_rows[2], y_rows[3], ty);
    }

    return CatmullRom(z_rows[0], z_rows[1], z_rows[2], z_rows[3], tz);
}

Real SampleTricubicLagrange(
    const std::vector<Real>& data,
    int width,
    int height,
    int depth,
    Real x,
    Real y,
    Real z) {
    const Real cx = std::clamp(x, 0.0, static_cast<Real>(width - 1));
    const Real cy = std::clamp(y, 0.0, static_cast<Real>(height - 1));
    const Real cz = std::clamp(z, 0.0, static_cast<Real>(depth - 1));
    const int ix = static_cast<int>(std::floor(cx));
    const int iy = static_cast<int>(std::floor(cy));
    const int iz = static_cast<int>(std::floor(cz));

    Real wx[4];
    Real wy[4];
    Real wz[4];
    LagrangeWeights(cx - static_cast<Real>(ix), wx);
    LagrangeWeights(cy - static_cast<Real>(iy), wy);
    LagrangeWeights(cz - static_cast<Real>(iz), wz);

    Real result = 0.0;
    for (int kz = 0; kz < 4; ++kz) {
        for (int jy = 0; jy < 4; ++jy) {
            for (int ix_offset = 0; ix_offset < 4; ++ix_offset) {
                result += wx[ix_offset] * wy[jy] * wz[kz] *
                          GetAt(data, width, height, depth, ix + ix_offset - 1, iy + jy - 1, iz + kz - 1);
            }
        }
    }

    return result;
}

}  // namespace skyspaces::reference
