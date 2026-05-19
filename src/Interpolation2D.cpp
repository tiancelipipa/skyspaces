#include "Interpolation2D.h"

#include <algorithm>
#include <cmath>

namespace skyspaces {

namespace {

// Clamp helper: constrain `v` to the closed interval [a, b]
inline double Clamp(double v, double a, double b) {
    return std::max(a, std::min(v, b));
}

// Sample a scalar grid stored in row-major order (width x height).
inline double GetAt(const std::vector<double>& d, int w, int h, int x, int y) {
    x = std::max(0, std::min(x, w - 1));
    y = std::max(0, std::min(y, h - 1));
    return d[static_cast<std::size_t>(y) * static_cast<std::size_t>(w) + static_cast<std::size_t>(x)];
}

// Catmull-Rom cubic interpolation (1D)
inline double CatmullRom(double p0, double p1, double p2, double p3, double t) {
    double a = -0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3;
    double b = p0 - 2.5*p1 + 2.0*p2 - 0.5*p3;
    double c = -0.5*p0 + 0.5*p2;
    double d = p1;
    return ((a*t + b)*t + c)*t + d;
}

// Lagrange cubic weight for four-point interpolation at parameter s in [0,1].
// Uses points q_{i-1}, q_i, q_{i+1}, q_{i+2}.
inline void LagrangeWeights(double s, double weights[4]) {
    double s2 = s * s;
    double s3 = s2 * s;
    weights[0] = -s3 / 6.0 + s2 / 2.0 - s / 3.0;
    weights[1] = 1.0 - 0.5 * s - s2 + 0.5 * s3;
    weights[2] = s + 0.5 * s2 - 0.5 * s3;
    weights[3] = s3 / 6.0 - s / 6.0;
}

} // anonymous

double SampleNearest(const std::vector<double>& data, int width, int height, double x, double y) {
    // if (data.empty() || width <= 0 || height <= 0) return 0.0;
    int xi = static_cast<int>(std::round(Clamp(x, 0.0, static_cast<double>(width - 1))));
    int yi = static_cast<int>(std::round(Clamp(y, 0.0, static_cast<double>(height - 1))));
    return GetAt(data, width, height, xi, yi);
}

double SampleBilinear(const std::vector<double>& data, int width, int height, double x, double y) {
    // if (data.empty() || width <= 0 || height <= 0) return 0.0;
    double cx = Clamp(x, 0.0, static_cast<double>(width - 1));
    double cy = Clamp(y, 0.0, static_cast<double>(height - 1));
    int x0 = static_cast<int>(std::floor(cx));
    int y0 = static_cast<int>(std::floor(cy));
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    double tx = cx - x0;
    double ty = cy - y0;

    double c00 = GetAt(data, width, height, x0, y0);
    double c10 = GetAt(data, width, height, x1, y0);
    double c01 = GetAt(data, width, height, x0, y1);
    double c11 = GetAt(data, width, height, x1, y1);

    double c0 = (1.0 - tx) * c00 + tx * c10;
    double c1 = (1.0 - tx) * c01 + tx * c11;
    return (1.0 - ty) * c0 + ty * c1;
}

double SampleBicubicCatmullRom(const std::vector<double>& data, int width, int height, double x, double y) {
    // if (data.empty() || width <= 0 || height <= 0) return 0.0;
    double cx = Clamp(x, 0.0, static_cast<double>(width - 1));
    double cy = Clamp(y, 0.0, static_cast<double>(height - 1));
    int ix = static_cast<int>(std::floor(cx));
    int iy = static_cast<int>(std::floor(cy));
    double tx = cx - ix;
    double ty = cy - iy;

    // Gather 4x4 neighborhood
    double cols[4];
    double tmp[4];
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

double SampleBicubicLagrange(const std::vector<double>& data, int width, int height, double x, double y) {
    // if (data.empty() || width <= 0 || height <= 0) return 0.0;
    double cx = Clamp(x, 0.0, static_cast<double>(width - 1));
    double cy = Clamp(y, 0.0, static_cast<double>(height - 1));
    int ix = static_cast<int>(std::floor(cx));
    int iy = static_cast<int>(std::floor(cy));
    double tx = cx - ix;
    double ty = cy - iy;

    double wx[4];
    double wy[4];
    LagrangeWeights(tx, wx);
    LagrangeWeights(ty, wy);

    double row[4];
    for (int j = 0; j < 4; ++j) {
        int yy = iy + j - 1;
        double v0 = GetAt(data, width, height, ix - 1, yy);
        double v1 = GetAt(data, width, height, ix, yy);
        double v2 = GetAt(data, width, height, ix + 1, yy);
        double v3 = GetAt(data, width, height, ix + 2, yy);
        row[j] = wx[0] * v0 + wx[1] * v1 + wx[2] * v2 + wx[3] * v3;
    }

    return wy[0] * row[0] + wy[1] * row[1] + wy[2] * row[2] + wy[3] * row[3];
}

} // namespace skyspaces
