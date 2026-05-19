#include "MACGrid2D.h"
#include "Interpolation2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace skyspaces {

MACScalarGrid2D::MACScalarGrid2D(int width, int height, double value) {
    Resize(width, height, value);
}

void MACScalarGrid2D::Resize(int width, int height, double value) {
    assert(width >= 0 && height >= 0);
    width_ = width;
    height_ = height;
    data_.assign(static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_), value);
}

void MACScalarGrid2D::Fill(double value) {
    std::fill(data_.begin(), data_.end(), value);
}

int MACScalarGrid2D::Width() const noexcept {
    return width_;
}

int MACScalarGrid2D::Height() const noexcept {
    return height_;
}

bool MACScalarGrid2D::Empty() const noexcept {
    return data_.empty();
}

double& MACScalarGrid2D::operator()(int x, int y) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    return data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) + static_cast<std::size_t>(x)];
}

double MACScalarGrid2D::operator()(int x, int y) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    return data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) + static_cast<std::size_t>(x)];
}

double MACScalarGrid2D::Sample(double x, double y) const {
    return SampleBilinear(data_, width_, height_, x, y);
}

const std::vector<double>& MACScalarGrid2D::Data() const noexcept {
    return data_;
}

std::vector<double>& MACScalarGrid2D::Data() noexcept {
    return data_;
}

MACVectorGrid2D::MACVectorGrid2D(int width, int height, double value) {
    Resize(width, height, value);
}

void MACVectorGrid2D::Resize(int width, int height, double value) {
    assert(width >= 0 && height >= 0);
    width_ = width;
    height_ = height;
    u_data_.assign(static_cast<std::size_t>(width_ + 1) * static_cast<std::size_t>(height_), value);
    v_data_.assign(static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_ + 1), value);
}

void MACVectorGrid2D::Fill(double value) {
    std::fill(u_data_.begin(), u_data_.end(), value);
    std::fill(v_data_.begin(), v_data_.end(), value);
}

int MACVectorGrid2D::Width() const noexcept {
    return width_;
}

int MACVectorGrid2D::Height() const noexcept {
    return height_;
}

int MACVectorGrid2D::UWidth() const noexcept {
    return width_ + 1;
}

int MACVectorGrid2D::UHeight() const noexcept {
    return height_;
}

int MACVectorGrid2D::VWidth() const noexcept {
    return width_;
}

int MACVectorGrid2D::VHeight() const noexcept {
    return height_ + 1;
}

bool MACVectorGrid2D::Empty() const noexcept {
    return u_data_.empty() && v_data_.empty();
}

double& MACVectorGrid2D::operator()(int x, int y, int component) {
    assert(component == 0 || component == 1);
    if (component == 0) {
        return U(x, y);
    }
    return V(x, y);
}

double MACVectorGrid2D::operator()(int x, int y, int component) const {
    assert(component == 0 || component == 1);
    if (component == 0) {
        return U(x, y);
    }
    return V(x, y);
}

double& MACVectorGrid2D::U(int x, int y) {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    return u_data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_ + 1) + static_cast<std::size_t>(x)];
}

double MACVectorGrid2D::U(int x, int y) const {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    return u_data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_ + 1) + static_cast<std::size_t>(x)];
}

double& MACVectorGrid2D::V(int x, int y) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    return v_data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) + static_cast<std::size_t>(x)];
}

double MACVectorGrid2D::V(int x, int y) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    return v_data_[static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) + static_cast<std::size_t>(x)];
}

Vector2D MACVectorGrid2D::Sample(double x, double y) const {
    const double u = SampleBilinear(u_data_, width_ + 1, height_, x, y);  
    const double v = SampleBilinear(v_data_, width_, height_ + 1, x, y);  
    return Vector2D(u, v);
}

const std::vector<double>& MACVectorGrid2D::UData() const noexcept {
    return u_data_;
}

std::vector<double>& MACVectorGrid2D::UData() noexcept {
    return u_data_;
}

const std::vector<double>& MACVectorGrid2D::VData() const noexcept {
    return v_data_;
}

std::vector<double>& MACVectorGrid2D::VData() noexcept {
    return v_data_;
}

}  // namespace skyspaces

