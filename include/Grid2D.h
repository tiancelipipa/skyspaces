#pragma once

#include "Interpolation2D.h"
#include "MathTypes.h"

namespace skyspaces {

// Scalar quantities such as smoke density, temperature, pressure, divergence,
// and vorticity live at cell centers. Public indexing is (x, y), and the
// Eigen array also stores rows as x and columns as y.
class CellCenteredScalarGrid2D {
public:
    CellCenteredScalarGrid2D() = default;
    CellCenteredScalarGrid2D(int width, int height, Real value = 0.0);

    void Resize(int width, int height, Real value = 0.0);
    void Fill(Real value);

    int Width() const noexcept;
    int Height() const noexcept;
    bool Empty() const noexcept;

    Real& operator()(int x, int y);
    Real operator()(int x, int y) const;

    // Samples in grid-index coordinates, where cell centers are integer points.
    Real Sample(Real x, Real y) const;
    Real Sample(Real x, Real y, InterpolationMethod2D method) const;
    Real Sample(const Vector2D& position) const;
    Real Sample(const Vector2D& position, InterpolationMethod2D method) const;

    const ScalarArray2D& Data() const noexcept;
    ScalarArray2D& Data() noexcept;

private:
    int width_ = 0;
    int height_ = 0;
    ScalarArray2D data_;
};

// Staggered/MAC velocity grid. Horizontal velocity U is stored on vertical
// faces with size (width + 1, height); vertical velocity V is stored on
// horizontal faces with size (width, height + 1).
class FaceCenteredVectorGrid2D {
public:
    FaceCenteredVectorGrid2D() = default;
    FaceCenteredVectorGrid2D(int width, int height, Real value = 0.0);

    void Resize(int width, int height, Real value = 0.0);
    void Fill(Real value);

    int Width() const noexcept;
    int Height() const noexcept;
    int UWidth() const noexcept;
    int UHeight() const noexcept;
    int VWidth() const noexcept;
    int VHeight() const noexcept;
    bool Empty() const noexcept;

    Real& operator()(int x, int y, int component);
    Real operator()(int x, int y, int component) const;

    Real& U(int x, int y);
    Real U(int x, int y) const;
    Real& V(int x, int y);
    Real V(int x, int y) const;

    // Samples both components in their own staggered index spaces.
    Vector2D Sample(Real x, Real y) const;
    Vector2D Sample(Real x, Real y, InterpolationMethod2D method) const;
    Vector2D Sample(const Vector2D& position) const;
    Vector2D Sample(const Vector2D& position, InterpolationMethod2D method) const;
    Real SampleU(Real x, Real y) const;
    Real SampleU(Real x, Real y, InterpolationMethod2D method) const;
    Real SampleV(Real x, Real y) const;
    Real SampleV(Real x, Real y, InterpolationMethod2D method) const;

    const ScalarArray2D& UData() const noexcept;
    ScalarArray2D& UData() noexcept;
    const ScalarArray2D& VData() const noexcept;
    ScalarArray2D& VData() noexcept;

private:
    int width_ = 0;
    int height_ = 0;
    ScalarArray2D u_data_;
    ScalarArray2D v_data_;
};

// Transitional aliases for older call sites. Prefer the explicit names above
// in new code, especially if a 3D grid is added later.
using MACScalarGrid2D = CellCenteredScalarGrid2D;
using MACVectorGrid2D = FaceCenteredVectorGrid2D;

}  // namespace skyspaces
