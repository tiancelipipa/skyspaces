#include "Solid2D.h"

#include "GridUtils2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace skyspaces {

namespace {

const CellCenteredScalarGrid2D& EmptyScalarGrid2D() {
    static const CellCenteredScalarGrid2D grid;
    return grid;
}

const FaceCenteredVectorGrid2D& EmptyVectorGrid2D() {
    static const FaceCenteredVectorGrid2D grid;
    return grid;
}

}  // namespace

Solid2D::Solid2D(int resolution_x, int resolution_y, Real cell_size, Real numeric_epsilon) {
    Resize(resolution_x, resolution_y, cell_size, numeric_epsilon);
}

void Solid2D::Resize(int resolution_x, int resolution_y, Real cell_size, Real numeric_epsilon) {
    assert(resolution_x > 0 && resolution_y > 0);
    assert(cell_size > 0.0);
    assert(numeric_epsilon > 0.0);

    resolution_x_ = resolution_x;
    resolution_y_ = resolution_y;
    cell_size_ = cell_size;
    inverse_cell_size_ = 1.0 / cell_size_;
    numeric_epsilon_ = numeric_epsilon;

    cell_markers_.reset();
    level_set_.reset();
    velocity_.reset();
    mode_ = SolidBoundaryMode2D::None;
}

void Solid2D::Clear() {
    mode_ = SolidBoundaryMode2D::None;
    cell_markers_.reset();
    level_set_.reset();
    velocity_.reset();
}

void Solid2D::SetCellMarkers(const CellCenteredScalarGrid2D& solid_cells) {
    if (solid_cells.Width() != resolution_x_ || solid_cells.Height() != resolution_y_) {
        throw std::invalid_argument("solid cell marker grid resolution does not match Solid2D");
    }

    if (!cell_markers_) {
        cell_markers_.emplace(resolution_x_, resolution_y_, 0.0);
    }
    cell_markers_->Data() = solid_cells.Data();
    level_set_.reset();
    mode_ = SolidBoundaryMode2D::CellMarkers;
}

void Solid2D::SetLevelSet(const CellCenteredScalarGrid2D& solid_phi) {
    if (solid_phi.Width() != resolution_x_ || solid_phi.Height() != resolution_y_) {
        throw std::invalid_argument("solid level set grid resolution does not match Solid2D");
    }

    if (!level_set_) {
        level_set_.emplace(resolution_x_, resolution_y_, 1.0);
    }
    level_set_->Data() = solid_phi.Data();
    cell_markers_.reset();
    mode_ = SolidBoundaryMode2D::LevelSet;
}

void Solid2D::SetVelocity(const FaceCenteredVectorGrid2D& solid_velocity) {
    if (solid_velocity.Width() != resolution_x_ || solid_velocity.Height() != resolution_y_) {
        throw std::invalid_argument("solid velocity grid resolution does not match Solid2D");
    }

    if (!velocity_) {
        velocity_.emplace(resolution_x_, resolution_y_, 0.0);
    }
    velocity_->UData() = solid_velocity.UData();
    velocity_->VData() = solid_velocity.VData();
}

int Solid2D::ResolutionX() const noexcept {
    return resolution_x_;
}

int Solid2D::ResolutionY() const noexcept {
    return resolution_y_;
}

Real Solid2D::CellSize() const noexcept {
    return cell_size_;
}

SolidBoundaryMode2D Solid2D::Mode() const noexcept {
    return mode_;
}

bool Solid2D::HasBoundary() const noexcept {
    return mode_ != SolidBoundaryMode2D::None;
}

bool Solid2D::IsCellSolid(int x, int y) const {
    if (!HasBoundary()) {
        return false;
    }

    if (x < 0 || x >= resolution_x_ || y < 0 || y >= resolution_y_) {
        return false;
    }

    switch (mode_) {
        case SolidBoundaryMode2D::CellMarkers:
            return cell_markers_ && (*cell_markers_)(x, y) > 0.5;
        case SolidBoundaryMode2D::LevelSet:
            return level_set_ && (*level_set_)(x, y) < 0.0;
        case SolidBoundaryMode2D::None:
            return false;
    }

    return false;
}

bool Solid2D::IsWorldSolid(Real x, Real y) const {
    if (!HasBoundary()) {
        return false;
    }

    if (mode_ == SolidBoundaryMode2D::LevelSet) {
        return SampleLevelSetWorld(x, y) < 0.0;
    }

    const int cell_x = ClampGridIndex2D(
        static_cast<int>(std::floor(x * inverse_cell_size_)),
        resolution_x_);
    const int cell_y = ClampGridIndex2D(
        static_cast<int>(std::floor(y * inverse_cell_size_)),
        resolution_y_);
    return IsCellSolid(cell_x, cell_y);
}

bool Solid2D::IsUFaceOpen(int x, int y) const {
    if (!HasBoundary()) {
        return true;
    }

    const bool left_is_fluid = x > 0 && !IsCellSolid(x - 1, y);
    const bool right_is_fluid = x < resolution_x_ && !IsCellSolid(x, y);
    return left_is_fluid && right_is_fluid;
}

bool Solid2D::IsVFaceOpen(int x, int y) const {
    if (!HasBoundary()) {
        return true;
    }

    const bool bottom_is_fluid = y > 0 && !IsCellSolid(x, y - 1);
    const bool top_is_fluid = y < resolution_y_ && !IsCellSolid(x, y);
    return bottom_is_fluid && top_is_fluid;
}

bool Solid2D::IsUFaceAdjacentToSolid(int x, int y) const {
    if (!HasBoundary()) {
        return false;
    }

    const bool left_is_solid = x > 0 && IsCellSolid(x - 1, y);
    const bool right_is_solid = x < resolution_x_ && IsCellSolid(x, y);
    return left_is_solid || right_is_solid;
}

bool Solid2D::IsVFaceAdjacentToSolid(int x, int y) const {
    if (!HasBoundary()) {
        return false;
    }

    const bool bottom_is_solid = y > 0 && IsCellSolid(x, y - 1);
    const bool top_is_solid = y < resolution_y_ && IsCellSolid(x, y);
    return bottom_is_solid || top_is_solid;
}

Real Solid2D::UFaceVelocity(int x, int y) const {
    if (!velocity_ || x < 0 || x >= velocity_->UWidth() || y < 0 || y >= velocity_->UHeight()) {
        return 0.0;
    }
    return IsUFaceAdjacentToSolid(x, y) ? velocity_->U(x, y) : 0.0;
}

Real Solid2D::VFaceVelocity(int x, int y) const {
    if (!velocity_ || x < 0 || x >= velocity_->VWidth() || y < 0 || y >= velocity_->VHeight()) {
        return 0.0;
    }
    return IsVFaceAdjacentToSolid(x, y) ? velocity_->V(x, y) : 0.0;
}

Vector2R Solid2D::ProjectOutOfSolid(const Vector2R& position, const Vector2R& fallback) const {
    if (!IsWorldSolid(position.x(), position.y())) {
        return position;
    }

    if (!IsWorldSolid(fallback.x(), fallback.y())) {
        if (mode_ == SolidBoundaryMode2D::LevelSet) {
            return ProjectOutOfLevelSetSolid_(position, fallback);
        }
        return ProjectOutOfMarkedSolid_(position, fallback);
    }

    return position;
}

Real Solid2D::SampleLevelSetWorld(Real x, Real y) const {
    if (!level_set_) {
        return 1.0;
    }

    const Real grid_x = CellCenterGridCoordinate2D(x, inverse_cell_size_);
    const Real grid_y = CellCenterGridCoordinate2D(y, inverse_cell_size_);
    return level_set_->Sample(grid_x, grid_y, InterpolationMethod2D::Bilinear);
}

Vector2R Solid2D::LevelSetGradientWorld(Real x, Real y) const {
    const Real half_dx = 0.5 * cell_size_;
    const Real max_x = static_cast<Real>(resolution_x_) * cell_size_;
    const Real max_y = static_cast<Real>(resolution_y_) * cell_size_;
    const Real x_left = std::clamp(x - half_dx, 0.0, max_x);
    const Real x_right = std::clamp(x + half_dx, 0.0, max_x);
    const Real y_bottom = std::clamp(y - half_dx, 0.0, max_y);
    const Real y_top = std::clamp(y + half_dx, 0.0, max_y);

    const Real dx = std::max(x_right - x_left, numeric_epsilon_);
    const Real dy = std::max(y_top - y_bottom, numeric_epsilon_);
    return {
        (SampleLevelSetWorld(x_right, y) - SampleLevelSetWorld(x_left, y)) / dx,
        (SampleLevelSetWorld(x, y_top) - SampleLevelSetWorld(x, y_bottom)) / dy};
}

const CellCenteredScalarGrid2D& Solid2D::CellMarkers() const noexcept {
    return cell_markers_ ? *cell_markers_ : EmptyScalarGrid2D();
}

const CellCenteredScalarGrid2D& Solid2D::LevelSet() const noexcept {
    return level_set_ ? *level_set_ : EmptyScalarGrid2D();
}

const FaceCenteredVectorGrid2D& Solid2D::Velocity() const noexcept {
    return velocity_ ? *velocity_ : EmptyVectorGrid2D();
}

Vector2R Solid2D::ProjectOutOfLevelSetSolid_(
    const Vector2R& position,
    const Vector2R& fallback) const {
    const Vector2R domain_max{
        static_cast<Real>(resolution_x_) * cell_size_,
        static_cast<Real>(resolution_y_) * cell_size_};
    Vector2R projected = position;

    for (int iter = 0; iter < 8; ++iter) {
        const Real phi = SampleLevelSetWorld(projected.x(), projected.y());
        if (phi >= 0.0) {
            break;
        }

        Vector2R normal = LevelSetGradientWorld(projected.x(), projected.y());
        const Real normal_length = normal.norm();
        if (normal_length <= numeric_epsilon_) {
            return fallback;
        }
        normal /= normal_length;
        projected -= phi * normal;
        projected.x() = std::clamp(projected.x(), 0.0, domain_max.x());
        projected.y() = std::clamp(projected.y(), 0.0, domain_max.y());
    }

    return IsWorldSolid(projected.x(), projected.y()) ? fallback : projected;
}

Vector2R Solid2D::ProjectOutOfMarkedSolid_(
    const Vector2R& position,
    const Vector2R& fallback) const {
    const Vector2R domain_max{
        static_cast<Real>(resolution_x_) * cell_size_,
        static_cast<Real>(resolution_y_) * cell_size_};
    Vector2R projected = position;
    Vector2R direction = fallback - position;
    if (direction.squaredNorm() <= numeric_epsilon_) {
        direction = {cell_size_, 0.0};
    }

    direction.normalize();
    const Real max_distance = (fallback - position).norm() + 2.0 * cell_size_;
    const int steps = std::max(1, static_cast<int>(std::ceil(max_distance / cell_size_ * 2.0)));
    const Real step = max_distance / static_cast<Real>(steps);
    for (int i = 0; i <= steps; ++i) {
        projected = position + direction * (static_cast<Real>(i) * step);
        projected.x() = std::clamp(projected.x(), 0.0, domain_max.x());
        projected.y() = std::clamp(projected.y(), 0.0, domain_max.y());
        if (!IsWorldSolid(projected.x(), projected.y())) {
            return projected;
        }
    }

    return fallback;
}

}  // namespace skyspaces
