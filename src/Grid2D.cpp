#include "Grid2D.h"

#include "Interpolation2D.h"

#include <cassert>

namespace skyspaces {

CellCenteredScalarGrid2D::CellCenteredScalarGrid2D(int width, int height, Real value) {
    Resize(width, height, value);
}

void CellCenteredScalarGrid2D::Resize(int width, int height, Real value) {
    assert(width >= 0 && height >= 0);
    width_ = width;
    height_ = height;
    data_.resize(width_, height_);
    data_.setConstant(value);
}

void CellCenteredScalarGrid2D::Fill(Real value) {
    data_.setConstant(value);
}

int CellCenteredScalarGrid2D::Width() const noexcept {
    return width_;
}

int CellCenteredScalarGrid2D::Height() const noexcept {
    return height_;
}

bool CellCenteredScalarGrid2D::Empty() const noexcept {
    return data_.size() == 0;
}

Real& CellCenteredScalarGrid2D::operator()(int x, int y) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    return data_(x, y);
}

Real CellCenteredScalarGrid2D::operator()(int x, int y) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    return data_(x, y);
}

Real CellCenteredScalarGrid2D::Sample(Real x, Real y) const {
    return SampleBilinear(data_, width_, height_, x, y);
}

Real CellCenteredScalarGrid2D::Sample(Real x, Real y, InterpolationMethod2D method) const {
    return skyspaces::Sample(data_, width_, height_, x, y, method);
}

Real CellCenteredScalarGrid2D::Sample(const Vector2R& position) const {
    return Sample(position.x(), position.y());
}

Real CellCenteredScalarGrid2D::Sample(const Vector2R& position, InterpolationMethod2D method) const {
    return Sample(position.x(), position.y(), method);
}

const ScalarArray2D& CellCenteredScalarGrid2D::Data() const noexcept {
    return data_;
}

ScalarArray2D& CellCenteredScalarGrid2D::Data() noexcept {
    return data_;
}

FaceCenteredVectorGrid2D::FaceCenteredVectorGrid2D(int width, int height, Real value) {
    Resize(width, height, value);
}

void FaceCenteredVectorGrid2D::Resize(int width, int height, Real value) {
    assert(width >= 0 && height >= 0);
    width_ = width;
    height_ = height;
    u_data_.resize(width_ + 1, height_);
    u_data_.setConstant(value);
    v_data_.resize(width_, height_ + 1);
    v_data_.setConstant(value);
}

void FaceCenteredVectorGrid2D::Fill(Real value) {
    u_data_.setConstant(value);
    v_data_.setConstant(value);
}

int FaceCenteredVectorGrid2D::Width() const noexcept {
    return width_;
}

int FaceCenteredVectorGrid2D::Height() const noexcept {
    return height_;
}

int FaceCenteredVectorGrid2D::UWidth() const noexcept {
    return width_ + 1;
}

int FaceCenteredVectorGrid2D::UHeight() const noexcept {
    return height_;
}

int FaceCenteredVectorGrid2D::VWidth() const noexcept {
    return width_;
}

int FaceCenteredVectorGrid2D::VHeight() const noexcept {
    return height_ + 1;
}

bool FaceCenteredVectorGrid2D::Empty() const noexcept {
    return u_data_.size() == 0 && v_data_.size() == 0;
}

Real& FaceCenteredVectorGrid2D::operator()(int x, int y, int component) {
    assert(component == 0 || component == 1);
    return component == 0 ? U(x, y) : V(x, y);
}

Real FaceCenteredVectorGrid2D::operator()(int x, int y, int component) const {
    assert(component == 0 || component == 1);
    return component == 0 ? U(x, y) : V(x, y);
}

Real& FaceCenteredVectorGrid2D::U(int x, int y) {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    return u_data_(x, y);
}

Real FaceCenteredVectorGrid2D::U(int x, int y) const {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    return u_data_(x, y);
}

Real& FaceCenteredVectorGrid2D::V(int x, int y) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    return v_data_(x, y);
}

Real FaceCenteredVectorGrid2D::V(int x, int y) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    return v_data_(x, y);
}

Vector2R FaceCenteredVectorGrid2D::Sample(Real x, Real y) const {
    return {SampleU(x, y), SampleV(x, y)};
}

Vector2R FaceCenteredVectorGrid2D::Sample(Real x, Real y, InterpolationMethod2D method) const {
    return {SampleU(x, y, method), SampleV(x, y, method)};
}

Vector2R FaceCenteredVectorGrid2D::Sample(const Vector2R& position) const {
    return Sample(position.x(), position.y());
}

Vector2R FaceCenteredVectorGrid2D::Sample(const Vector2R& position, InterpolationMethod2D method) const {
    return Sample(position.x(), position.y(), method);
}

Real FaceCenteredVectorGrid2D::SampleU(Real x, Real y) const {
    return SampleBilinear(u_data_, width_ + 1, height_, x, y);
}

Real FaceCenteredVectorGrid2D::SampleU(Real x, Real y, InterpolationMethod2D method) const {
    return skyspaces::Sample(u_data_, width_ + 1, height_, x, y, method);
}

Real FaceCenteredVectorGrid2D::SampleV(Real x, Real y) const {
    return SampleBilinear(v_data_, width_, height_ + 1, x, y);
}

Real FaceCenteredVectorGrid2D::SampleV(Real x, Real y, InterpolationMethod2D method) const {
    return skyspaces::Sample(v_data_, width_, height_ + 1, x, y, method);
}

const ScalarArray2D& FaceCenteredVectorGrid2D::UData() const noexcept {
    return u_data_;
}

ScalarArray2D& FaceCenteredVectorGrid2D::UData() noexcept {
    return u_data_;
}

const ScalarArray2D& FaceCenteredVectorGrid2D::VData() const noexcept {
    return v_data_;
}

ScalarArray2D& FaceCenteredVectorGrid2D::VData() noexcept {
    return v_data_;
}

}  // namespace skyspaces
