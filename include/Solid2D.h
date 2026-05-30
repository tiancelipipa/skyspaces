#pragma once

#include "Grid2D.h"

#include <optional>

namespace skyspaces {

enum class SolidBoundaryMode2D {
    None,
    CellMarkers,
    LevelSet,
};

// Stores the 2D solid boundary representation used by the fluid solver.
// Geometry can be supplied either as cell markers (> 0.5 means solid) or as a
// cell-centered signed distance field (negative means solid). Moving solids use
// a staggered velocity grid sampled on the same faces as the fluid velocity.
class Solid2D {
public:
    Solid2D() = default;
    Solid2D(int resolution_x, int resolution_y, Real cell_size, Real numeric_epsilon = 1e-12);

    void Resize(int resolution_x, int resolution_y, Real cell_size, Real numeric_epsilon = 1e-12);
    void Clear();

    void SetCellMarkers(const CellCenteredScalarGrid2D& solid_cells);
    void SetLevelSet(const CellCenteredScalarGrid2D& solid_phi);
    void SetVelocity(const FaceCenteredVectorGrid2D& solid_velocity);

    int ResolutionX() const noexcept;
    int ResolutionY() const noexcept;
    Real CellSize() const noexcept;
    SolidBoundaryMode2D Mode() const noexcept;
    bool HasBoundary() const noexcept;

    bool IsCellSolid(int x, int y) const;
    bool IsWorldSolid(Real x, Real y) const;
    bool IsUFaceOpen(int x, int y) const;
    bool IsVFaceOpen(int x, int y) const;
    bool IsUFaceAdjacentToSolid(int x, int y) const;
    bool IsVFaceAdjacentToSolid(int x, int y) const;
    Real UFaceVelocity(int x, int y) const;
    Real VFaceVelocity(int x, int y) const;

    Vector2R ProjectOutOfSolid(const Vector2R& position, const Vector2R& fallback) const;
    Real SampleLevelSetWorld(Real x, Real y) const;
    Vector2R LevelSetGradientWorld(Real x, Real y) const;

    const CellCenteredScalarGrid2D& CellMarkers() const noexcept;
    const CellCenteredScalarGrid2D& LevelSet() const noexcept;
    const FaceCenteredVectorGrid2D& Velocity() const noexcept;

private:
    Vector2R ProjectOutOfLevelSetSolid_(const Vector2R& position, const Vector2R& fallback) const;
    Vector2R ProjectOutOfMarkedSolid_(const Vector2R& position, const Vector2R& fallback) const;

    int resolution_x_ = 0;
    int resolution_y_ = 0;
    Real cell_size_ = 1.0;
    Real inverse_cell_size_ = 1.0;
    Real numeric_epsilon_ = 1e-12;
    SolidBoundaryMode2D mode_ = SolidBoundaryMode2D::None;

    std::optional<CellCenteredScalarGrid2D> cell_markers_;
    std::optional<CellCenteredScalarGrid2D> level_set_;
    std::optional<FaceCenteredVectorGrid2D> velocity_;
};

}  // namespace skyspaces
