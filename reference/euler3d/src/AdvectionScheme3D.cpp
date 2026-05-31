#include "AdvectionScheme3D.h"

#include "EulerFluid3D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <utility>

namespace skyspaces::reference {
namespace {

Real CellCenterGridCoordinate(Real world, Real inverse_cell_size) {
    return world * inverse_cell_size - 0.5;
}

Real FaceGridCoordinate(Real world, Real inverse_cell_size) {
    return world * inverse_cell_size;
}

template <typename ValueAt>
std::pair<Real, Real> GridValueRange3D(
    int width,
    int height,
    int depth,
    Real grid_x,
    Real grid_y,
    Real grid_z,
    ValueAt value_at) {
    const int x0 = static_cast<int>(std::floor(grid_x));
    const int y0 = static_cast<int>(std::floor(grid_y));
    const int z0 = static_cast<int>(std::floor(grid_z));
    Real min_value = std::numeric_limits<Real>::infinity();
    Real max_value = -std::numeric_limits<Real>::infinity();

    for (int dx = 0; dx <= 1; ++dx) {
        for (int dy = 0; dy <= 1; ++dy) {
            for (int dz = 0; dz <= 1; ++dz) {
                const int x = std::clamp(x0 + dx, 0, width - 1);
                const int y = std::clamp(y0 + dy, 0, height - 1);
                const int z = std::clamp(z0 + dz, 0, depth - 1);
                const Real value = value_at(x, y, z);
                min_value = std::min(min_value, value);
                max_value = std::max(max_value, value);
            }
        }
    }

    return {min_value, max_value};
}

}  // namespace

struct EulerFluid3DAdvection3D {
    static void AdvectVelocitySemiLagrangianScheme(EulerFluid3D& fluid, Real dt) {
        AdvectVelocitySemiLagrangian(fluid, fluid.velocity_, fluid.velocity_tmp_, dt);
    }

    static void AdvectVelocityMacCormackBFECCScheme(EulerFluid3D& fluid, Real dt) {
        assert(fluid.velocity_first_pass_);
        assert(fluid.velocity_back_pass_);
        assert(fluid.velocity_corrected_source_);
        FaceCenteredVectorGrid3D& first_pass = *fluid.velocity_first_pass_;
        FaceCenteredVectorGrid3D& back_pass = *fluid.velocity_back_pass_;
        FaceCenteredVectorGrid3D& corrected_source = *fluid.velocity_corrected_source_;

        AdvectVelocitySemiLagrangian(fluid, fluid.velocity_, first_pass, dt);
        AdvectVelocitySemiLagrangian(fluid, first_pass, back_pass, -dt);

        for (int x = 0; x < fluid.velocity_.UWidth(); ++x) {
            for (int y = 0; y < fluid.velocity_.UHeight(); ++y) {
                for (int z = 0; z < fluid.velocity_.UDepth(); ++z) {
                    const Real world_x = static_cast<Real>(x) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const auto range = GridValueRange3D(
                        fluid.velocity_.UWidth(),
                        fluid.velocity_.UHeight(),
                        fluid.velocity_.UDepth(),
                        FaceGridCoordinate(world_x, fluid.inverse_cell_size_),
                        CellCenterGridCoordinate(world_y, fluid.inverse_cell_size_),
                        CellCenterGridCoordinate(world_z, fluid.inverse_cell_size_),
                        [&](int ix, int iy, int iz) {
                            return fluid.velocity_.U(ix, iy, iz);
                        });
                    corrected_source.U(x, y, z) = std::clamp(
                        fluid.velocity_.U(x, y, z) +
                            0.5 * (fluid.velocity_.U(x, y, z) - back_pass.U(x, y, z)),
                        range.first,
                        range.second);
                }
            }
        }
        for (int x = 0; x < fluid.velocity_.VWidth(); ++x) {
            for (int y = 0; y < fluid.velocity_.VHeight(); ++y) {
                for (int z = 0; z < fluid.velocity_.VDepth(); ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = static_cast<Real>(y) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const auto range = GridValueRange3D(
                        fluid.velocity_.VWidth(),
                        fluid.velocity_.VHeight(),
                        fluid.velocity_.VDepth(),
                        CellCenterGridCoordinate(world_x, fluid.inverse_cell_size_),
                        FaceGridCoordinate(world_y, fluid.inverse_cell_size_),
                        CellCenterGridCoordinate(world_z, fluid.inverse_cell_size_),
                        [&](int ix, int iy, int iz) {
                            return fluid.velocity_.V(ix, iy, iz);
                        });
                    corrected_source.V(x, y, z) = std::clamp(
                        fluid.velocity_.V(x, y, z) +
                            0.5 * (fluid.velocity_.V(x, y, z) - back_pass.V(x, y, z)),
                        range.first,
                        range.second);
                }
            }
        }
        for (int x = 0; x < fluid.velocity_.WWidth(); ++x) {
            for (int y = 0; y < fluid.velocity_.WHeight(); ++y) {
                for (int z = 0; z < fluid.velocity_.WDepth(); ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = static_cast<Real>(z) * fluid.config_.cell_size;
                    const auto range = GridValueRange3D(
                        fluid.velocity_.WWidth(),
                        fluid.velocity_.WHeight(),
                        fluid.velocity_.WDepth(),
                        CellCenterGridCoordinate(world_x, fluid.inverse_cell_size_),
                        CellCenterGridCoordinate(world_y, fluid.inverse_cell_size_),
                        FaceGridCoordinate(world_z, fluid.inverse_cell_size_),
                        [&](int ix, int iy, int iz) {
                            return fluid.velocity_.W(ix, iy, iz);
                        });
                    corrected_source.W(x, y, z) = std::clamp(
                        fluid.velocity_.W(x, y, z) +
                            0.5 * (fluid.velocity_.W(x, y, z) - back_pass.W(x, y, z)),
                        range.first,
                        range.second);
                }
            }
        }

        AdvectVelocitySemiLagrangian(fluid, corrected_source, fluid.velocity_tmp_, dt);
    }

    static void AdvectScalarsSemiLagrangianScheme(EulerFluid3D& fluid, Real dt) {
        AdvectScalarSemiLagrangian(fluid, fluid.smoke_density_, fluid.smoke_density_tmp_, dt);
        AdvectScalarSemiLagrangian(fluid, fluid.temperature_, fluid.temperature_tmp_, dt);
    }

    static void AdvectScalarsMacCormackBFECCScheme(EulerFluid3D& fluid, Real dt) {
        assert(fluid.smoke_first_pass_);
        assert(fluid.temperature_first_pass_);
        assert(fluid.smoke_back_pass_);
        assert(fluid.temperature_back_pass_);
        assert(fluid.smoke_corrected_source_);
        assert(fluid.temperature_corrected_source_);
        CellCenteredScalarGrid3D& smoke_first = *fluid.smoke_first_pass_;
        CellCenteredScalarGrid3D& temperature_first = *fluid.temperature_first_pass_;
        CellCenteredScalarGrid3D& smoke_back = *fluid.smoke_back_pass_;
        CellCenteredScalarGrid3D& temperature_back = *fluid.temperature_back_pass_;
        CellCenteredScalarGrid3D& smoke_corrected = *fluid.smoke_corrected_source_;
        CellCenteredScalarGrid3D& temperature_corrected = *fluid.temperature_corrected_source_;

        AdvectScalarSemiLagrangian(fluid, fluid.smoke_density_, smoke_first, dt);
        AdvectScalarSemiLagrangian(fluid, fluid.temperature_, temperature_first, dt);
        AdvectScalarSemiLagrangian(fluid, smoke_first, smoke_back, -dt);
        AdvectScalarSemiLagrangian(fluid, temperature_first, temperature_back, -dt);

        for (int x = 0; x < fluid.config_.resolution_x; ++x) {
            for (int y = 0; y < fluid.config_.resolution_y; ++y) {
                for (int z = 0; z < fluid.config_.resolution_z; ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const Real grid_x = CellCenterGridCoordinate(world_x, fluid.inverse_cell_size_);
                    const Real grid_y = CellCenterGridCoordinate(world_y, fluid.inverse_cell_size_);
                    const Real grid_z = CellCenterGridCoordinate(world_z, fluid.inverse_cell_size_);
                    const auto smoke_range = GridValueRange3D(
                        fluid.smoke_density_.Width(),
                        fluid.smoke_density_.Height(),
                        fluid.smoke_density_.Depth(),
                        grid_x,
                        grid_y,
                        grid_z,
                        [&](int ix, int iy, int iz) {
                            return fluid.smoke_density_(ix, iy, iz);
                        });
                    const auto temperature_range = GridValueRange3D(
                        fluid.temperature_.Width(),
                        fluid.temperature_.Height(),
                        fluid.temperature_.Depth(),
                        grid_x,
                        grid_y,
                        grid_z,
                        [&](int ix, int iy, int iz) {
                            return fluid.temperature_(ix, iy, iz);
                        });
                    smoke_corrected(x, y, z) = std::clamp(
                        fluid.smoke_density_(x, y, z) +
                            0.5 * (fluid.smoke_density_(x, y, z) - smoke_back(x, y, z)),
                        smoke_range.first,
                        smoke_range.second);
                    temperature_corrected(x, y, z) = std::clamp(
                        fluid.temperature_(x, y, z) +
                            0.5 * (fluid.temperature_(x, y, z) - temperature_back(x, y, z)),
                        temperature_range.first,
                        temperature_range.second);
                }
            }
        }

        AdvectScalarSemiLagrangian(fluid, smoke_corrected, fluid.smoke_density_tmp_, dt);
        AdvectScalarSemiLagrangian(fluid, temperature_corrected, fluid.temperature_tmp_, dt);
    }

private:
    static void AdvectVelocitySemiLagrangian(
        const EulerFluid3D& fluid,
        const FaceCenteredVectorGrid3D& source,
        FaceCenteredVectorGrid3D& target,
        Real dt) {
        for (int x = 0; x < target.UWidth(); ++x) {
            for (int y = 0; y < target.UHeight(); ++y) {
                for (int z = 0; z < target.UDepth(); ++z) {
                    const Real world_x = static_cast<Real>(x) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const Vector3R sample_position = fluid.TraceAdvectionPosition_({world_x, world_y, world_z}, dt);
                    target.U(x, y, z) =
                        fluid.SampleUWorld_(source, sample_position.x(), sample_position.y(), sample_position.z());
                }
            }
        }

        for (int x = 0; x < target.VWidth(); ++x) {
            for (int y = 0; y < target.VHeight(); ++y) {
                for (int z = 0; z < target.VDepth(); ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = static_cast<Real>(y) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const Vector3R sample_position = fluid.TraceAdvectionPosition_({world_x, world_y, world_z}, dt);
                    target.V(x, y, z) =
                        fluid.SampleVWorld_(source, sample_position.x(), sample_position.y(), sample_position.z());
                }
            }
        }

        for (int x = 0; x < target.WWidth(); ++x) {
            for (int y = 0; y < target.WHeight(); ++y) {
                for (int z = 0; z < target.WDepth(); ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = static_cast<Real>(z) * fluid.config_.cell_size;
                    const Vector3R sample_position = fluid.TraceAdvectionPosition_({world_x, world_y, world_z}, dt);
                    target.W(x, y, z) =
                        fluid.SampleWWorld_(source, sample_position.x(), sample_position.y(), sample_position.z());
                }
            }
        }
    }

    static void AdvectScalarSemiLagrangian(
        const EulerFluid3D& fluid,
        const CellCenteredScalarGrid3D& source,
        CellCenteredScalarGrid3D& target,
        Real dt) {
        for (int x = 0; x < fluid.config_.resolution_x; ++x) {
            for (int y = 0; y < fluid.config_.resolution_y; ++y) {
                for (int z = 0; z < fluid.config_.resolution_z; ++z) {
                    const Real world_x = (static_cast<Real>(x) + 0.5) * fluid.config_.cell_size;
                    const Real world_y = (static_cast<Real>(y) + 0.5) * fluid.config_.cell_size;
                    const Real world_z = (static_cast<Real>(z) + 0.5) * fluid.config_.cell_size;
                    const Vector3R sample_position = fluid.TraceAdvectionPosition_({world_x, world_y, world_z}, dt);
                    target(x, y, z) = fluid.SampleScalarWorld_(
                        source,
                        sample_position.x(),
                        sample_position.y(),
                        sample_position.z());
                }
            }
        }
    }
};

void AdvectVelocitySemiLagrangian3D(EulerFluid3D& fluid, Real dt) {
    EulerFluid3DAdvection3D::AdvectVelocitySemiLagrangianScheme(fluid, dt);
}

void AdvectScalarsSemiLagrangian3D(EulerFluid3D& fluid, Real dt) {
    EulerFluid3DAdvection3D::AdvectScalarsSemiLagrangianScheme(fluid, dt);
}

void AdvectVelocityMacCormackBFECC3D(EulerFluid3D& fluid, Real dt) {
    EulerFluid3DAdvection3D::AdvectVelocityMacCormackBFECCScheme(fluid, dt);
}

void AdvectScalarsMacCormackBFECC3D(EulerFluid3D& fluid, Real dt) {
    EulerFluid3DAdvection3D::AdvectScalarsMacCormackBFECCScheme(fluid, dt);
}

}  // namespace skyspaces::reference
