#pragma once

#include "Vector2D.h"
#include <vector>

namespace skyspaces {

class MACScalarGrid2D {
public:
    MACScalarGrid2D() = default;
    MACScalarGrid2D(int width, int height, double value = 0.0);

    void Resize(int width, int height, double value = 0.0);
    void Fill(double value);

    int Width() const noexcept;
    int Height() const noexcept;
    bool Empty() const noexcept;

    double& operator()(int x, int y);
    double operator()(int x, int y) const;

    double Sample(double x, double y) const;

    const std::vector<double>& Data() const noexcept;
    std::vector<double>& Data() noexcept;

private:
    int width_ = 0;     // number of cells in x direction
    int height_ = 0;    // number of cells in y direction
    std::vector<double> data_;  // size: width_ * height_
};

class MACVectorGrid2D {
public:
    MACVectorGrid2D() = default;
    MACVectorGrid2D(int width, int height, double value = 0.0);

    void Resize(int width, int height, double value = 0.0);
    void Fill(double value);

    int Width() const noexcept;
    int Height() const noexcept;
    int UWidth() const noexcept;
    int UHeight() const noexcept;
    int VWidth() const noexcept;
    int VHeight() const noexcept;
    bool Empty() const noexcept;

    double& operator()(int x, int y, int component);
    double operator()(int x, int y, int component) const;

    double& U(int x, int y);
    double U(int x, int y) const;
    double& V(int x, int y);
    double V(int x, int y) const;

    Vector2D Sample(double x, double y) const;

    const std::vector<double>& UData() const noexcept;
    std::vector<double>& UData() noexcept;
    const std::vector<double>& VData() const noexcept;
    std::vector<double>& VData() noexcept;

private:
    int width_ = 0;     // number of cells in x direction
    int height_ = 0;    // number of cells in y direction
    std::vector<double> u_data_; // size: (width_ + 1) * height_
    std::vector<double> v_data_; // size: width_ * (height_ + 1)
};

}  // namespace skyspaces
