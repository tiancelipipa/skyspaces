#include "AdvectionScheme2D.h"

#include "AdvectionIntegrator2D.h"
#include "Fluid2D.h"
#include "GridUtils2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

namespace skyspaces {
namespace {

template <typename ValueAt>
std::pair<Real, Real> GridValueRange2D(
    int width,
    int height,
    Real grid_x,
    Real grid_y,
    ValueAt value_at) {
    const int x0 = static_cast<int>(std::floor(grid_x));
    const int y0 = static_cast<int>(std::floor(grid_y));
    Real min_value = std::numeric_limits<Real>::infinity();
    Real max_value = -std::numeric_limits<Real>::infinity();

    for (int dx = 0; dx <= 1; ++dx) {
        for (int dy = 0; dy <= 1; ++dy) {
            const int x = ClampGridIndex2D(x0 + dx, width);
            const int y = ClampGridIndex2D(y0 + dy, height);
            const Real value = value_at(x, y);
            min_value = std::min(min_value, value);
            max_value = std::max(max_value, value);
        }
    }

    return {min_value, max_value};
}

std::pair<Real, Real> CellCenteredRangeWorld2D(
    const CellCenteredScalarGrid2D& grid,
    Real world_x,
    Real world_y,
    Real inverse_cell_size) {
    return GridValueRange2D(
        grid.Width(),
        grid.Height(),
        CellCenterGridCoordinate2D(world_x, inverse_cell_size),
        CellCenterGridCoordinate2D(world_y, inverse_cell_size),
        [&](int x, int y) {
            return grid(x, y);
        });
}

std::pair<Real, Real> UFaceRangeWorld2D(
    const FaceCenteredVectorGrid2D& grid,
    Real world_x,
    Real world_y,
    Real inverse_cell_size) {
    return GridValueRange2D(
        grid.UWidth(),
        grid.UHeight(),
        FaceGridCoordinate2D(world_x, inverse_cell_size),
        CellCenterGridCoordinate2D(world_y, inverse_cell_size),
        [&](int x, int y) {
            return grid.U(x, y);
        });
}

std::pair<Real, Real> VFaceRangeWorld2D(
    const FaceCenteredVectorGrid2D& grid,
    Real world_x,
    Real world_y,
    Real inverse_cell_size) {
    return GridValueRange2D(
        grid.VWidth(),
        grid.VHeight(),
        CellCenterGridCoordinate2D(world_x, inverse_cell_size),
        FaceGridCoordinate2D(world_y, inverse_cell_size),
        [&](int x, int y) {
            return grid.V(x, y);
        });
}

}  // namespace

struct Fluid2DAdvection2D {
    static void AdvectVelocitySemiLagrangianScheme(Fluid2D& fluid, Real dt) {
        AdvectVelocitySemiLagrangian(fluid, fluid.velocity_, fluid.velocity_tmp_, dt);
    }

    static void AdvectVelocityMacCormackBFECCScheme(Fluid2D& fluid, Real dt) {
        assert(fluid.velocity_first_pass_);
        assert(fluid.velocity_back_pass_);
        assert(fluid.velocity_corrected_source_);
        FaceCenteredVectorGrid2D& first_pass = *fluid.velocity_first_pass_;
        FaceCenteredVectorGrid2D& back_pass = *fluid.velocity_back_pass_;
        FaceCenteredVectorGrid2D& corrected_source = *fluid.velocity_corrected_source_;

        AdvectVelocitySemiLagrangian(fluid, fluid.velocity_, first_pass, dt);
        AdvectVelocitySemiLagrangian(fluid, first_pass, back_pass, -dt);

        for (int x = 0; x < fluid.velocity_.UWidth(); ++x) {
            for (int y = 0; y < fluid.velocity_.UHeight(); ++y) {
                const Real world_x = static_cast<Real>(x) * fluid.config_.cell_size;
                const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                const auto range =
                    UFaceRangeWorld2D(fluid.velocity_, world_x, world_y, fluid.inverse_cell_size_);
                corrected_source.U(x, y) = std::clamp(
                    fluid.velocity_.U(x, y) + 0.5 * (fluid.velocity_.U(x, y) - back_pass.U(x, y)),
                    range.first,
                    range.second);
            }
        }
        for (int x = 0; x < fluid.velocity_.VWidth(); ++x) {
            for (int y = 0; y < fluid.velocity_.VHeight(); ++y) {
                const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                const Real world_y = static_cast<Real>(y) * fluid.config_.cell_size;
                const auto range =
                    VFaceRangeWorld2D(fluid.velocity_, world_x, world_y, fluid.inverse_cell_size_);
                corrected_source.V(x, y) = std::clamp(
                    fluid.velocity_.V(x, y) + 0.5 * (fluid.velocity_.V(x, y) - back_pass.V(x, y)),
                    range.first,
                    range.second);
            }
        }
        AdvectVelocitySemiLagrangian(fluid, corrected_source, fluid.velocity_tmp_, dt);
    }

    static void AdvectScalarsSemiLagrangianScheme(Fluid2D& fluid, Real dt) {
        AdvectScalarSemiLagrangian(fluid, fluid.smoke_density_, fluid.smoke_density_tmp_, dt, 0.0);
        AdvectScalarSemiLagrangian(
            fluid,
            fluid.temperature_,
            fluid.temperature_tmp_,
            dt,
            fluid.config_.ambient_temperature);
    }

    static void AdvectScalarsMacCormackBFECCScheme(Fluid2D& fluid, Real dt) {
        assert(fluid.smoke_first_pass_);
        assert(fluid.temperature_first_pass_);
        assert(fluid.smoke_back_pass_);
        assert(fluid.temperature_back_pass_);
        assert(fluid.smoke_corrected_source_);
        assert(fluid.temperature_corrected_source_);
        CellCenteredScalarGrid2D& smoke_first = *fluid.smoke_first_pass_;
        CellCenteredScalarGrid2D& temperature_first = *fluid.temperature_first_pass_;
        CellCenteredScalarGrid2D& smoke_back = *fluid.smoke_back_pass_;
        CellCenteredScalarGrid2D& temperature_back = *fluid.temperature_back_pass_;
        CellCenteredScalarGrid2D& smoke_corrected = *fluid.smoke_corrected_source_;
        CellCenteredScalarGrid2D& temperature_corrected = *fluid.temperature_corrected_source_;

        AdvectScalarSemiLagrangian(fluid, fluid.smoke_density_, smoke_first, dt, 0.0);
        AdvectScalarSemiLagrangian(
            fluid,
            fluid.temperature_,
            temperature_first,
            dt,
            fluid.config_.ambient_temperature);
        AdvectScalarSemiLagrangian(fluid, smoke_first, smoke_back, -dt, 0.0);
        AdvectScalarSemiLagrangian(
            fluid,
            temperature_first,
            temperature_back,
            -dt,
            fluid.config_.ambient_temperature);

        for (int x = 0; x < fluid.config_.resolution_x; ++x) {
            for (int y = 0; y < fluid.config_.resolution_y; ++y) {
                const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                const auto smoke_range =
                    CellCenteredRangeWorld2D(fluid.smoke_density_, world_x, world_y, fluid.inverse_cell_size_);
                const auto temperature_range =
                    CellCenteredRangeWorld2D(fluid.temperature_, world_x, world_y, fluid.inverse_cell_size_);
                smoke_corrected(x, y) = std::clamp(
                    fluid.smoke_density_(x, y) +
                        0.5 * (fluid.smoke_density_(x, y) - smoke_back(x, y)),
                    smoke_range.first,
                    smoke_range.second);
                temperature_corrected(x, y) = std::clamp(
                    fluid.temperature_(x, y) +
                        0.5 * (fluid.temperature_(x, y) - temperature_back(x, y)),
                    temperature_range.first,
                    temperature_range.second);
            }
        }

        AdvectScalarSemiLagrangian(fluid, smoke_corrected, fluid.smoke_density_tmp_, dt, 0.0);
        AdvectScalarSemiLagrangian(
            fluid,
            temperature_corrected,
            fluid.temperature_tmp_,
            dt,
            fluid.config_.ambient_temperature);
    }

private:
    static Vector2R TraceAdvectionPosition(const Fluid2D& fluid, const Vector2R& position, Real dt) {
        const auto velocity_at = [&](const Vector2R& p) {
            return fluid.SampleVelocityWorld_(fluid.velocity_, p.x(), p.y());
        };
        return fluid.ProjectOutOfSolid_(
            BacktracePosition2D(fluid.config_.advection_integrator, position, dt, velocity_at),
            position);
    }

    static void AdvectVelocitySemiLagrangian(
        const Fluid2D& fluid,
        const FaceCenteredVectorGrid2D& source,
        FaceCenteredVectorGrid2D& target,
        Real dt) {
        for (int x = 0; x < target.UWidth(); ++x) {
            for (int y = 0; y < target.UHeight(); ++y) {
                const Real world_x = static_cast<Real>(x) * fluid.config_.cell_size;
                const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                const Vector2R sample_position = TraceAdvectionPosition(fluid, {world_x, world_y}, dt);
                target.U(x, y) = fluid.IsUFaceOpen_(x, y)
                    ? fluid.SampleUWorld_(source, sample_position.x(), sample_position.y())
                    : fluid.SolidU_(x, y);
            }
        }

        for (int x = 0; x < target.VWidth(); ++x) {
            for (int y = 0; y < target.VHeight(); ++y) {
                const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                const Real world_y = static_cast<Real>(y) * fluid.config_.cell_size;
                const Vector2R sample_position = TraceAdvectionPosition(fluid, {world_x, world_y}, dt);
                target.V(x, y) = fluid.IsVFaceOpen_(x, y)
                    ? fluid.SampleVWorld_(source, sample_position.x(), sample_position.y())
                    : fluid.SolidV_(x, y);
            }
        }
    }

    static void AdvectScalarSemiLagrangian(
        const Fluid2D& fluid,
        const CellCenteredScalarGrid2D& source,
        CellCenteredScalarGrid2D& target,
        Real dt,
        Real solid_value) {
        for (int x = 0; x < fluid.config_.resolution_x; ++x) {
            for (int y = 0; y < fluid.config_.resolution_y; ++y) {
                if (fluid.IsSolidCell_(x, y)) {
                    target(x, y) = solid_value;
                    continue;
                }

                const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                const Vector2R sample_position = TraceAdvectionPosition(fluid, {world_x, world_y}, dt);
                target(x, y) = fluid.SampleScalarWorld_(source, sample_position.x(), sample_position.y());
            }
        }
    }
};

void AdvectVelocitySemiLagrangian2D(Fluid2D& fluid, Real dt) {
    Fluid2DAdvection2D::AdvectVelocitySemiLagrangianScheme(fluid, dt);
}

void AdvectScalarsSemiLagrangian2D(Fluid2D& fluid, Real dt) {
    Fluid2DAdvection2D::AdvectScalarsSemiLagrangianScheme(fluid, dt);
}

void AdvectVelocityMacCormackBFECC2D(Fluid2D& fluid, Real dt) {
    Fluid2DAdvection2D::AdvectVelocityMacCormackBFECCScheme(fluid, dt);
}

void AdvectScalarsMacCormackBFECC2D(Fluid2D& fluid, Real dt) {
    Fluid2DAdvection2D::AdvectScalarsMacCormackBFECCScheme(fluid, dt);
}

}  // namespace skyspaces
