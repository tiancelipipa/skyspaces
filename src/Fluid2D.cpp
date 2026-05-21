#include "Fluid2D.h"

#include <Eigen/Sparse>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace skyspaces {

namespace {

Real Clamp(Real value, Real min_value, Real max_value) {
    return std::max(min_value, std::min(value, max_value));
}

int ClampIndex(int value, int min_value, int max_value) {
    return std::max(min_value, std::min(value, max_value));
}

int CellIndex(int x, int y, int height) {
    // Keep vectorized pressure storage aligned with the row-major x/y grid.
    return x * height + y;
}

}  // namespace

Fluid2D::Fluid2D(Config config) : config_(config) {
    InitializeFromConfig_();
}

void Fluid2D::Reset() {
    smoke_density_.Fill(0.0);
    smoke_density_tmp_.Fill(0.0);
    temperature_.Fill(config_.ambient_temperature);
    temperature_tmp_.Fill(config_.ambient_temperature);
    pressure_.Fill(0.0);
    divergence_.Fill(0.0);
    vorticity_.Fill(0.0);
    velocity_.Fill(0.0);
    velocity_tmp_.Fill(0.0);

    time_ = 0.0;
    step_count_ = 0;
    last_pressure_iterations_ = 0;
    last_pressure_residual_ = 0.0;
}

void Fluid2D::Step() {
    Step(config_.time_step);
}

void Fluid2D::Step(Real dt) {
    if (dt <= 0.0) {
        return;
    }

    const int substeps = ComputeSubstepCount_(dt);
    const Real substep_dt = dt / static_cast<Real>(substeps);

    for (int substep = 0; substep < substeps; ++substep) {
        // Advect with the divergence-free velocity from the previous substep.
        AdvectVelocity_(substep_dt);

        // Body(External) forces and sources are applied before projection 
        // so that they affect the velocity field used for advection. 
        ApplyConfiguredVelocitySource_(substep_dt);
        ApplyBuoyancy_(substep_dt);

        if (config_.vorticity_confinement > 0.0) {
            ComputeVorticity_();
            ApplyVorticityConfinement_(substep_dt);
        }

        // Pressure projection makes the velocity field divergence-free(incompressible)
        SetBoundaryConditions_();
        ProjectVelocity_(substep_dt);

        // Advect smoke and temperature with the updated divergence-free velocity field.
        AdvectScalars_(substep_dt);
        ApplyConfiguredScalarSource_(substep_dt);

        time_ += substep_dt;
    }

    ++step_count_;
}

void Fluid2D::AddSmokeSourceNormalized_(
    Real min_x,
    Real max_x,
    Real min_y,
    Real max_y,
    Real smoke_rate,
    Real target_temperature,
    Real temperature_rate,
    Real dt) {

    assert (min_x >= 0.0 && min_x < max_x && max_x <= 1.0);
    assert (min_y >= 0.0 && min_y < max_y && max_y <= 1.0);
    assert (temperature_rate >= 0.0);
    
    const int x_begin = static_cast<int>(std::floor(min_x * config_.resolution_x));
    const int x_end = static_cast<int>(std::ceil(max_x * config_.resolution_x));

    const int y_begin = static_cast<int>(std::floor(min_y * config_.resolution_y));
    const int y_end = static_cast<int>(std::ceil(max_y * config_.resolution_y));

    const int width = x_end - x_begin;
    const int height = y_end - y_begin;

    auto smoke = smoke_density_.Data().block(x_begin, y_begin, width, height);
    smoke = (smoke + smoke_rate * dt).min(config_.max_smoke_density).max(0.0);
    
    const Real heat_alpha = 1.0 - std::exp(-std::max(0.0, temperature_rate) * dt);
    auto temperature = temperature_.Data().block(x_begin, y_begin, width, height);
    temperature += heat_alpha * (target_temperature - temperature);
}

void Fluid2D::AddVelocityImpulseNormalized_(
    Real min_x,
    Real max_x,
    Real min_y,
    Real max_y,
    Vector2D acceleration,
    Real dt) {
    assert (min_x >= 0.0 && min_x < max_x && max_x <= 1.0);
    assert (min_y >= 0.0 && min_y < max_y && max_y <= 1.0);
    
    const Real width = static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real height = static_cast<Real>(config_.resolution_y) * config_.cell_size;
    const Real velocity_delta_x = acceleration.x() * dt * width;
    const Real velocity_delta_y = acceleration.y() * dt * height;

    for (int x = 0; x < velocity_.UWidth(); ++x) {
        for (int y = 0; y < velocity_.UHeight(); ++y) {
            const Real py = (static_cast<Real>(y) + 0.5)  / static_cast<Real>(config_.resolution_y);
            if (py < min_y || py > max_y) {
                continue;
            }

            const Real px = static_cast<Real>(x) / static_cast<Real>(config_.resolution_x);
            if (px >= min_x && px <= max_x) {
                velocity_.U(x, y) += velocity_delta_x;
            }
        }
    }

    for (int x = 0; x < velocity_.VWidth(); ++x) {
        for (int y = 0; y < velocity_.VHeight(); ++y) {
            const Real py = static_cast<Real>(y) / static_cast<Real>(config_.resolution_y);
            if (py < min_y || py > max_y) {
                continue;
            }

            const Real px = (static_cast<Real>(x) + 0.5) / static_cast<Real>(config_.resolution_x);
            if (px >= min_x && px <= max_x) {
                velocity_.V(x, y) += velocity_delta_y;
            }
        }
    }
}

int Fluid2D::ResolutionX() const noexcept {
    return config_.resolution_x;
}

int Fluid2D::ResolutionY() const noexcept {
    return config_.resolution_y;
}

Real Fluid2D::CellSize() const noexcept {
    return config_.cell_size;
}

Real Fluid2D::TimeStep() const noexcept {
    return config_.time_step;
}

Real Fluid2D::Time() const noexcept {
    return time_;
}

int Fluid2D::StepCount() const noexcept {
    return step_count_;
}

int Fluid2D::LastPressureIterations() const noexcept {
    return last_pressure_iterations_;
}

Real Fluid2D::LastPressureResidual() const noexcept {
    return last_pressure_residual_;
}

Real Fluid2D::SampleSmokeNormalized(Real x, Real y) const {
    const Real world_x = Clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = Clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    return SampleScalarWorld_(smoke_density_, world_x, world_y);
}

Real Fluid2D::SampleTemperatureNormalized(Real x, Real y) const {
    const Real world_x = Clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = Clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    return SampleScalarWorld_(temperature_, world_x, world_y);
}

Vector2D Fluid2D::SampleVelocityNormalized(Real x, Real y) const {
    const Real world_x = Clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = Clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    return SampleVelocityWorld_(velocity_, world_x, world_y);
}

const Fluid2D::Config& Fluid2D::Configs() const noexcept {
    return config_;
}

Fluid2D::Config& Fluid2D::Configs() noexcept {
    return config_;
}

const CellCenteredScalarGrid2D& Fluid2D::SmokeDensity() const noexcept {
    return smoke_density_;
}

const CellCenteredScalarGrid2D& Fluid2D::Temperature() const noexcept {
    return temperature_;
}

const CellCenteredScalarGrid2D& Fluid2D::Pressure() const noexcept {
    return pressure_;
}

const CellCenteredScalarGrid2D& Fluid2D::Divergence() const noexcept {
    return divergence_;
}

const CellCenteredScalarGrid2D& Fluid2D::Vorticity() const noexcept {
    return vorticity_;
}

const FaceCenteredVectorGrid2D& Fluid2D::Velocity() const noexcept {
    return velocity_;
}

void Fluid2D::InitializeFromConfig_() {
    assert (config_.numeric_epsilon > 0.0);
    
    assert (config_.resolution_x > 0 && config_.resolution_y > 0);
    assert (config_.cell_size > 0.0);

    assert (config_.time_step > 0.0);
    assert (config_.max_substeps > 0);
    assert (config_.cfl_number >= 1.0 && config_.cfl_number <= 10.0);
    
    assert (config_.fluid_density > 0.0);
    assert (config_.pressure_iterations > 0);
    assert (config_.pressure_tolerance > 0.0);

    assert (config_.buoyancy_smoke_density_coefficient > 0.0);
    assert (config_.buoyancy_temperature_coefficient > 0.0);

    assert (config_.max_smoke_density > 0.0);
    
    assert (config_.source_min_x >= 0.0 && 
            config_.source_min_x < config_.source_max_x && 
            config_.source_max_x <= 1.0);
    assert (config_.source_min_y >= 0.0 && 
            config_.source_min_y < config_.source_max_y && 
            config_.source_max_y <= 1.0);
    assert (config_.source_temperature_rate >= 0.0);

    inverse_cell_size_ = 1.0 / config_.cell_size;

    smoke_density_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    smoke_density_tmp_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    temperature_.Resize(config_.resolution_x, config_.resolution_y, config_.ambient_temperature);
    temperature_tmp_.Resize(config_.resolution_x, config_.resolution_y, config_.ambient_temperature);
    pressure_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    divergence_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    vorticity_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    velocity_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
    velocity_tmp_.Resize(config_.resolution_x, config_.resolution_y, 0.0);
}

int Fluid2D::ComputeSubstepCount_(Real dt) const {
    const Real max_speed = ComputeMaxSpeed_();
    if (max_speed <= config_.numeric_epsilon) {
        return 1;
    }

    const Real max_dt = config_.cfl_number * config_.cell_size / max_speed;
    const int count = static_cast<int>(std::ceil(dt / std::max(max_dt, config_.numeric_epsilon)));
    return std::clamp(count, 1, config_.max_substeps);
}

Real Fluid2D::ComputeMaxSpeed_() const {
    const Real max_u = velocity_.UData().abs().maxCoeff();
    const Real max_v = velocity_.VData().abs().maxCoeff();
    return std::max(max_u, max_v);
}

void Fluid2D::ApplyConfiguredScalarSource_(Real dt) {
    if (!config_.source_enabled) {
        return;
    }

    AddSmokeSourceNormalized_(
        config_.source_min_x,
        config_.source_max_x,
        config_.source_min_y,
        config_.source_max_y,
        config_.source_smoke_rate,
        config_.source_temperature,
        config_.source_temperature_rate,
        dt);
}

void Fluid2D::ApplyConfiguredVelocitySource_(Real dt) {
    if (!config_.source_enabled) {
        return;
    }
    
    AddVelocityImpulseNormalized_(
        config_.source_min_x,
        config_.source_max_x,
        config_.source_min_y,
        config_.source_max_y,
        {0.0, config_.source_acceleration_y},
        dt);
}

void Fluid2D::ApplyBuoyancy_(Real dt) {
    const int width = config_.resolution_x;
    const int interior_height = config_.resolution_y - 1;
    if (interior_height <= 0) {
        return;
    }

    // V samples sit between scalar rows; average neighboring cell-centered
    // smoke and temperature before applying vertical acceleration.
    const ScalarArray2D smoke =
        0.5 * (smoke_density_.Data().block(0, 0, width, interior_height) +
               smoke_density_.Data().block(0, 1, width, interior_height)).eval();
    const ScalarArray2D temperature =
        0.5 * (temperature_.Data().block(0, 0, width, interior_height) +
               temperature_.Data().block(0, 1, width, interior_height)).eval();
    velocity_.VData().block(0, 1, width, interior_height) +=
        (config_.gravity *
         (config_.buoyancy_temperature_coefficient *
              (temperature - config_.ambient_temperature) -
          config_.buoyancy_smoke_density_coefficient * smoke)) *
        dt;
}

void Fluid2D::ComputeVorticity_() {
    const Real inv_dx = inverse_cell_size_;
    for (int x = 0; x < config_.resolution_x; ++x) {
        for (int y = 0; y < config_.resolution_y; ++y) {
            const Real world_x = (static_cast<Real>(x) + 0.5) * config_.cell_size;
            const Real world_y = (static_cast<Real>(y) + 0.5) * config_.cell_size;

            const Real dv_dx =
                (SampleVWorld_(velocity_, world_x + 0.5 * config_.cell_size, world_y) -
                 SampleVWorld_(velocity_, world_x - 0.5 * config_.cell_size, world_y)) *
                inv_dx;
            const Real du_dy =
                (SampleUWorld_(velocity_, world_x, world_y + 0.5 * config_.cell_size) -
                 SampleUWorld_(velocity_, world_x, world_y - 0.5 * config_.cell_size)) *
                inv_dx;
            vorticity_(x, y) = dv_dx - du_dy;
        }
    }
}

void Fluid2D::ApplyVorticityConfinement_(Real dt) {
    const Real eps = config_.vorticity_confinement;
    for (int x = 1; x < config_.resolution_x - 1; ++x) {
        for (int y = 1; y < config_.resolution_y - 1; ++y) {
            const Real grad_x = 0.5 * (std::abs(vorticity_(x + 1, y)) - std::abs(vorticity_(x - 1, y))) * inverse_cell_size_;
            const Real grad_y = 0.5 * (std::abs(vorticity_(x, y + 1)) - std::abs(vorticity_(x, y - 1))) * inverse_cell_size_;
            const Real len = std::sqrt(grad_x * grad_x + grad_y * grad_y);
            if (len <= config_.numeric_epsilon) {
                continue;
            }

            const Real nx = grad_x / len;
            const Real ny = grad_y / len;
            const Real curl = vorticity_(x, y);
            const Real force_x = eps * ny * curl;
            const Real force_y = -eps * nx * curl;

            velocity_.U(x, y) += force_x * dt;
            velocity_.U(x + 1, y) += force_x * dt;
            velocity_.V(x, y) += force_y * dt;
            velocity_.V(x, y + 1) += force_y * dt;
        }
    }
}

void Fluid2D::AdvectVelocity_(Real dt) {
    // Semi-Lagrangian advection: trace each face center backward through the
    // current velocity field, then sample the previous component there.
    const auto velocity_at = [&](const Vector2D& p) {
        return SampleVelocityWorld_(velocity_, p.x(), p.y());
    };

    for (int x = 0; x < velocity_.UWidth(); ++x) {
        for (int y = 0; y < velocity_.UHeight(); ++y) {
            const Real world_x = static_cast<Real>(x) * config_.cell_size;
            const Real world_y = (static_cast<Real>(y) + 0.5) * config_.cell_size;
            const Vector2D back_position =
                BacktracePosition2D(config_.advection_integrator, {world_x, world_y}, dt, velocity_at);
            velocity_tmp_.U(x, y) = SampleUWorld_(velocity_, back_position.x(), back_position.y());
        }
    }

    for (int x = 0; x < velocity_.VWidth(); ++x) {
        for (int y = 0; y < velocity_.VHeight(); ++y) {
            const Real world_x = (static_cast<Real>(x) + 0.5) * config_.cell_size;
            const Real world_y = static_cast<Real>(y) * config_.cell_size;
            const Vector2D back_position =
                BacktracePosition2D(config_.advection_integrator, {world_x, world_y}, dt, velocity_at);
            velocity_tmp_.V(x, y) = SampleVWorld_(velocity_, back_position.x(), back_position.y());
        }
    }

    velocity_.UData() = velocity_tmp_.UData();
    velocity_.VData() = velocity_tmp_.VData();
}

void Fluid2D::AdvectScalars_(Real dt) {
    const Real smoke_decay = std::exp(-std::max(0.0, config_.smoke_dissipation) * dt);
    const Real temp_decay = std::exp(-std::max(0.0, config_.temperature_dissipation) * dt);
    const auto velocity_at = [&](const Vector2D& p) {
        return SampleVelocityWorld_(velocity_, p.x(), p.y());
    };

    // Scalars are sampled at cell centers and decayed after advection.
    for (int x = 0; x < config_.resolution_x; ++x) {
        for (int y = 0; y < config_.resolution_y; ++y) {
            const Real world_x = (static_cast<Real>(x) + 0.5) * config_.cell_size;
            const Real world_y = (static_cast<Real>(y) + 0.5) * config_.cell_size;
            const Vector2D back_position =
                BacktracePosition2D(config_.advection_integrator, {world_x, world_y}, dt, velocity_at);

            smoke_density_tmp_(x, y) = Clamp(
                SampleScalarWorld_(smoke_density_, back_position.x(), back_position.y()) * smoke_decay,
                0.0,
                config_.max_smoke_density);

            const Real advected_temp =
                SampleScalarWorld_(temperature_, back_position.x(), back_position.y());
            temperature_tmp_(x, y) =
                config_.ambient_temperature +
                (advected_temp - config_.ambient_temperature) * temp_decay;
        }
    }

    smoke_density_.Data() = smoke_density_tmp_.Data();
    temperature_.Data() = temperature_tmp_.Data();
}

void Fluid2D::ProjectVelocity_(Real dt) {
    if (dt <= 0.0) {
        return;
    }

    ComputeDivergence_();

    const int width = config_.resolution_x;
    const int height = config_.resolution_y;
    const int cell_count = width * height;
    const Real scale = dt / (config_.fluid_density * config_.cell_size * config_.cell_size);

    // Build a 5-point Poisson stencil. Pinning cell 0 removes the constant
    // pressure nullspace.
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(static_cast<std::size_t>(cell_count) * 5);
    Eigen::VectorXd rhs(cell_count);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            const int row = CellIndex(x, y, height);
            if (row == 0) {
                // Fix one pressure value to define the pressure reference level.
                // With closed boundaries, pressure is only determined up to a
                // constant because velocity correction uses pressure gradients.
                triplets.emplace_back(row, row, 1.0);
                rhs(row) = 0.0;
                continue;
            }

            rhs(row) = -divergence_(x, y);

            Real diagonal = 0.0;
            const int offsets[4][2] = {
                {-1, 0},
                {1, 0},
                {0, -1},
                {0, 1},
            };

            for (const auto& offset : offsets) {
                const int nx = x + offset[0];
                const int ny = y + offset[1];
                if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
                    continue;
                }

                diagonal += scale;
                const int neighbor = CellIndex(nx, ny, height);
                if (neighbor != 0) {
                    triplets.emplace_back(row, neighbor, -scale);
                }
            }

            triplets.emplace_back(row, row, diagonal > 0.0 ? diagonal : 1.0);
        }
    }

    Eigen::SparseMatrix<Real> matrix(cell_count, cell_count);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();

    const PressureSolveResult2D pressure_result = SolvePressurePoisson2D(
        config_.pressure_solver,
        matrix,
        rhs,
        config_.pressure_iterations,
        config_.pressure_tolerance);

    if (!pressure_result.success) {
        last_pressure_iterations_ = 0;
        last_pressure_residual_ = std::numeric_limits<Real>::infinity();
        pressure_.Fill(0.0);
        return;
    }

    const Real max_pressure = pressure_result.solution.cwiseAbs().maxCoeff();
    const Real pressure_limit =
        std::max<Real>(1.0, rhs.cwiseAbs().maxCoeff()) * 1.0e6;
    if (!pressure_result.solution.allFinite() || max_pressure > pressure_limit) {
        last_pressure_iterations_ = pressure_result.iterations;
        last_pressure_residual_ = std::numeric_limits<Real>::infinity();
        pressure_.Fill(0.0);
        return;
    }

    last_pressure_iterations_ = pressure_result.iterations;
    last_pressure_residual_ = pressure_result.residual;

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            pressure_(x, y) = pressure_result.solution(CellIndex(x, y, height));
        }
    }

    const ScalarArray2D velocity_u_before_projection = velocity_.UData();
    const ScalarArray2D velocity_v_before_projection = velocity_.VData();
    const Real grad_scale = dt / (config_.fluid_density * config_.cell_size);
    // Subtract pressure gradients on interior faces, then clamp boundaries.
    velocity_.UData().block(1, 0, width - 1, height) -=
        grad_scale *
        (pressure_.Data().block(1, 0, width - 1, height) -
         pressure_.Data().block(0, 0, width - 1, height));
    velocity_.VData().block(0, 1, width, height - 1) -=
        grad_scale *
        (pressure_.Data().block(0, 1, width, height - 1) -
         pressure_.Data().block(0, 0, width, height - 1));

    SetBoundaryConditions_();
    ComputeDivergence_();

    const Real max_velocity_after_projection =
        std::max(velocity_.UData().abs().maxCoeff(), velocity_.VData().abs().maxCoeff());
    if (!std::isfinite(max_velocity_after_projection) || max_velocity_after_projection > 1.0e6) {
        velocity_.UData() = velocity_u_before_projection;
        velocity_.VData() = velocity_v_before_projection;
        pressure_.Fill(0.0);
        last_pressure_residual_ = std::numeric_limits<Real>::infinity();
        SetBoundaryConditions_();
        ComputeDivergence_();
    }
}

void Fluid2D::SetBoundaryConditions_() {
    velocity_.UData().row(0).setZero();
    velocity_.UData().row(velocity_.UWidth() - 1).setZero();
    velocity_.VData().col(0).setZero();
    velocity_.VData().col(velocity_.VHeight() - 1).setZero();
}

void Fluid2D::ComputeDivergence_() {
    const int width = config_.resolution_x;
    const int height = config_.resolution_y;
    divergence_.Data() =
        (velocity_.UData().block(1, 0, width, height) -
         velocity_.UData().block(0, 0, width, height) +
         velocity_.VData().block(0, 1, width, height) -
         velocity_.VData().block(0, 0, width, height)) *
        inverse_cell_size_;
}

Real Fluid2D::SampleScalarWorld_(const CellCenteredScalarGrid2D& grid, Real x, Real y) const {
    // World-space cell center (i + 0.5) * dx maps to grid index i.
    const Real grid_x = x * inverse_cell_size_ - 0.5;
    const Real grid_y = y * inverse_cell_size_ - 0.5;
    return grid.Sample(grid_x, grid_y, config_.advection_interpolation);
}

Vector2D Fluid2D::SampleVelocityWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    return {SampleUWorld_(grid, x, y), SampleVWorld_(grid, x, y)};
}

Real Fluid2D::SampleUWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    // U is face-aligned in x and cell-centered in y.
    const Real grid_x = x * inverse_cell_size_;
    const Real grid_y = y * inverse_cell_size_ - 0.5;
    return grid.SampleU(grid_x, grid_y, config_.advection_interpolation);
}

Real Fluid2D::SampleVWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    // V is cell-centered in x and face-aligned in y.
    const Real grid_x = x * inverse_cell_size_ - 0.5;
    const Real grid_y = y * inverse_cell_size_;
    return grid.SampleV(grid_x, grid_y, config_.advection_interpolation);
}

}  // namespace skyspaces
