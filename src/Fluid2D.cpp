#include "Fluid2D.h"

#include <Eigen/Sparse>

#include "GridUtils2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace skyspaces {

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
    solid_.Clear();
    ConfigureAdvectionScratch_();

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

    const Real heat_alpha = 1.0 - std::exp(-std::max(0.0, temperature_rate) * dt);
    for (int x = x_begin; x < x_begin + width; ++x) {
        for (int y = y_begin; y < y_begin + height; ++y) {
            if (IsSolidCell_(x, y)) {
                continue;
            }

            smoke_density_(x, y) = std::clamp(
                smoke_density_(x, y) + smoke_rate * dt,
                0.0,
                config_.max_smoke_density);
            temperature_(x, y) += heat_alpha * (target_temperature - temperature_(x, y));
        }
    }
}

void Fluid2D::AddVelocityImpulseNormalized_(
    Real min_x,
    Real max_x,
    Real min_y,
    Real max_y,
    Vector2R acceleration,
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
                if (!IsUFaceOpen_(x, y)) {
                    continue;
                }
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
                if (!IsVFaceOpen_(x, y)) {
                    continue;
                }
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

void Fluid2D::ClearSolidBoundary() {
    solid_.Clear();
    SetBoundaryConditions_();
}

void Fluid2D::SetSolidCellMarkers(const CellCenteredScalarGrid2D& solid_cells) {
    solid_.SetCellMarkers(solid_cells);
    SetBoundaryConditions_();
}

void Fluid2D::SetSolidLevelSet(const CellCenteredScalarGrid2D& solid_phi) {
    solid_.SetLevelSet(solid_phi);
    SetBoundaryConditions_();
}

void Fluid2D::SetSolidVelocity(const FaceCenteredVectorGrid2D& solid_velocity) {
    solid_.SetVelocity(solid_velocity);
    SetBoundaryConditions_();
}

Real Fluid2D::SampleSmokeNormalized(Real x, Real y) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    return SampleScalarWorld_(smoke_density_, world_x, world_y);
}

Real Fluid2D::SampleTemperatureNormalized(Real x, Real y) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    return SampleScalarWorld_(temperature_, world_x, world_y);
}

Vector2R Fluid2D::SampleVelocityNormalized(Real x, Real y) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
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

CellCenteredScalarGrid2D& Fluid2D::SmokeDensity() noexcept {
    return smoke_density_;
}

const CellCenteredScalarGrid2D& Fluid2D::Temperature() const noexcept {
    return temperature_;
}

CellCenteredScalarGrid2D& Fluid2D::Temperature() noexcept {
    return temperature_;
}

const CellCenteredScalarGrid2D& Fluid2D::Pressure() const noexcept {
    return pressure_;
}

CellCenteredScalarGrid2D& Fluid2D::Pressure() noexcept {
    return pressure_;
}

const CellCenteredScalarGrid2D& Fluid2D::Divergence() const noexcept {
    return divergence_;
}

const CellCenteredScalarGrid2D& Fluid2D::Vorticity() const noexcept {
    return vorticity_;
}

const CellCenteredScalarGrid2D& Fluid2D::SolidCellMarkers() const noexcept {
    return solid_.CellMarkers();
}

const CellCenteredScalarGrid2D& Fluid2D::SolidLevelSet() const noexcept {
    return solid_.LevelSet();
}

const FaceCenteredVectorGrid2D& Fluid2D::Velocity() const noexcept {
    return velocity_;
}

const FaceCenteredVectorGrid2D& Fluid2D::SolidVelocity() const noexcept {
    return solid_.Velocity();
}

const Solid2D& Fluid2D::Solid() const noexcept {
    return solid_;
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
    solid_.Resize(
        config_.resolution_x,
        config_.resolution_y,
        config_.cell_size,
        config_.numeric_epsilon);
    ConfigureAdvectionScratch_();
}

void Fluid2D::ConfigureAdvectionScratch_() {
    if (config_.advection_scheme == AdvectionScheme2D::SemiLagrangian) {
        ClearAdvectionScratch_();
        return;
    }

    velocity_first_pass_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    velocity_back_pass_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    velocity_corrected_source_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    smoke_first_pass_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    smoke_back_pass_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    smoke_corrected_source_.emplace(config_.resolution_x, config_.resolution_y, 0.0);
    temperature_first_pass_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.ambient_temperature);
    temperature_back_pass_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.ambient_temperature);
    temperature_corrected_source_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.ambient_temperature);
}

void Fluid2D::ClearAdvectionScratch_() {
    velocity_first_pass_.reset();
    velocity_back_pass_.reset();
    velocity_corrected_source_.reset();
    smoke_first_pass_.reset();
    smoke_back_pass_.reset();
    smoke_corrected_source_.reset();
    temperature_first_pass_.reset();
    temperature_back_pass_.reset();
    temperature_corrected_source_.reset();
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
        {config_.source_acceleration_x, config_.source_acceleration_y},
        dt);
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
    switch (config_.advection_scheme) {
    case AdvectionScheme2D::SemiLagrangian:
        AdvectVelocitySemiLagrangian2D(*this, dt);
        break;
    case AdvectionScheme2D::MacCormackBFECC:
        AdvectVelocityMacCormackBFECC2D(*this, dt);
        break;
    }

    velocity_.UData() = velocity_tmp_.UData();
    velocity_.VData() = velocity_tmp_.VData();
}

void Fluid2D::AdvectScalars_(Real dt) {
    switch (config_.advection_scheme) {
    case AdvectionScheme2D::SemiLagrangian:
        AdvectScalarsSemiLagrangian2D(*this, dt);
        break;
    case AdvectionScheme2D::MacCormackBFECC:
        AdvectScalarsMacCormackBFECC2D(*this, dt);
        break;
    }

    ApplyScalarAdvectionPostProcess_(dt);
}

void Fluid2D::ApplyScalarAdvectionPostProcess_(Real dt) {
    const Real smoke_decay = std::exp(-std::max(0.0, config_.smoke_dissipation) * dt);
    const Real temp_decay = std::exp(-std::max(0.0, config_.temperature_dissipation) * dt);

    for (int x = 0; x < config_.resolution_x; ++x) {
        for (int y = 0; y < config_.resolution_y; ++y) {
            if (IsSolidCell_(x, y)) {
                smoke_density_tmp_(x, y) = 0.0;
                temperature_tmp_(x, y) = config_.ambient_temperature;
                continue;
            }

            smoke_density_tmp_(x, y) = std::clamp(
                smoke_density_tmp_(x, y) * smoke_decay,
                0.0,
                config_.max_smoke_density);
            temperature_tmp_(x, y) =
                config_.ambient_temperature +
                (temperature_tmp_(x, y) - config_.ambient_temperature) * temp_decay;
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
    std::vector<int> component_id(static_cast<std::size_t>(cell_count), -1);
    std::vector<int> component_reference_rows;
    std::vector<Real> component_rhs_sums;
    std::vector<int> component_cell_counts;
    std::vector<int> stack;
    stack.reserve(static_cast<std::size_t>(cell_count));

    for (int start_x = 0; start_x < width; ++start_x) {
        for (int start_y = 0; start_y < height; ++start_y) {
            if (IsSolidCell_(start_x, start_y)) {
                continue;
            }

            const int start_row = FlattenCellIndex2D(start_x, start_y, height);
            if (component_id[static_cast<std::size_t>(start_row)] >= 0) {
                continue;
            }

            const int id = static_cast<int>(component_reference_rows.size());
            component_reference_rows.push_back(start_row);
            component_rhs_sums.push_back(0.0);
            component_cell_counts.push_back(0);

            stack.clear();
            stack.push_back(start_row);
            component_id[static_cast<std::size_t>(start_row)] = id;

            while (!stack.empty()) {
                const int row = stack.back();
                stack.pop_back();
                const int x = row / height;
                const int y = row - x * height;

                component_rhs_sums[static_cast<std::size_t>(id)] += -divergence_(x, y);
                ++component_cell_counts[static_cast<std::size_t>(id)];

                const int offsets[4][2] = {
                    {-1, 0},
                    {1, 0},
                    {0, -1},
                    {0, 1},
                };
                for (const auto& offset : offsets) {
                    const int nx = x + offset[0];
                    const int ny = y + offset[1];
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height || IsSolidCell_(nx, ny)) {
                        continue;
                    }

                    const int neighbor = FlattenCellIndex2D(nx, ny, height);
                    if (component_id[static_cast<std::size_t>(neighbor)] >= 0) {
                        continue;
                    }

                    component_id[static_cast<std::size_t>(neighbor)] = id;
                    stack.push_back(neighbor);
                }
            }
        }
    }

    if (component_reference_rows.empty()) {
        velocity_.Fill(0.0);
        pressure_.Fill(0.0);
        divergence_.Fill(0.0);
        return;
    }

    std::vector<char> is_pressure_reference(static_cast<std::size_t>(cell_count), 0);
    for (const int row : component_reference_rows) {
        is_pressure_reference[static_cast<std::size_t>(row)] = 1;
    }

    // Build a 5-point Poisson stencil. Pinning cell 0 removes the constant
    // pressure nullspace.
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(static_cast<std::size_t>(cell_count) * 5);
    Eigen::VectorXd rhs(cell_count);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            const int row = FlattenCellIndex2D(x, y, height);
            if (IsSolidCell_(x, y)) {
                triplets.emplace_back(row, row, 1.0);
                rhs(row) = 0.0;
                continue;
            }

            if (is_pressure_reference[static_cast<std::size_t>(row)]) {
                // Fix one pressure value to define the pressure reference level.
                // Each closed fluid component has its own constant-pressure
                // nullspace, so each component needs one pinned reference cell.
                triplets.emplace_back(row, row, 1.0);
                rhs(row) = 0.0;
                continue;
            }

            const int id = component_id[static_cast<std::size_t>(row)];
            const Real component_mean_rhs =
                component_rhs_sums[static_cast<std::size_t>(id)] /
                static_cast<Real>(component_cell_counts[static_cast<std::size_t>(id)]);
            rhs(row) = -divergence_(x, y) - component_mean_rhs;

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

                if (IsSolidCell_(nx, ny)) {
                    continue;
                }

                diagonal += scale;
                const int neighbor = FlattenCellIndex2D(nx, ny, height);
                if (!is_pressure_reference[static_cast<std::size_t>(neighbor)]) {
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
            pressure_(x, y) = IsSolidCell_(x, y)
                ? 0.0
                : pressure_result.solution(FlattenCellIndex2D(x, y, height));
        }
    }

    const ScalarArray2D velocity_u_before_projection = velocity_.UData();
    const ScalarArray2D velocity_v_before_projection = velocity_.VData();
    const Real grad_scale = dt / (config_.fluid_density * config_.cell_size);
    // Subtract pressure gradients on interior faces, then clamp boundaries.
    for (int x = 1; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (IsUFaceOpen_(x, y)) {
                velocity_.U(x, y) -= grad_scale * (pressure_(x, y) - pressure_(x - 1, y));
            } else {
                velocity_.U(x, y) = SolidU_(x, y);
            }
        }
    }
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            if (IsVFaceOpen_(x, y)) {
                velocity_.V(x, y) -= grad_scale * (pressure_(x, y) - pressure_(x, y - 1));
            } else {
                velocity_.V(x, y) = SolidV_(x, y);
            }
        }
    }

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

    if (!HasSolidBoundary_()) {
        return;
    }

    for (int x = 0; x < velocity_.UWidth(); ++x) {
        for (int y = 0; y < velocity_.UHeight(); ++y) {
            if (!IsUFaceOpen_(x, y)) {
                velocity_.U(x, y) = SolidU_(x, y);
            }
        }
    }

    for (int x = 0; x < velocity_.VWidth(); ++x) {
        for (int y = 0; y < velocity_.VHeight(); ++y) {
            if (!IsVFaceOpen_(x, y)) {
                velocity_.V(x, y) = SolidV_(x, y);
            }
        }
    }
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

    if (!HasSolidBoundary_()) {
        return;
    }

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            if (IsSolidCell_(x, y)) {
                divergence_(x, y) = 0.0;
            }
        }
    }
}

bool Fluid2D::HasSolidBoundary_() const noexcept {
    return solid_.HasBoundary();
}

bool Fluid2D::IsSolidCell_(int x, int y) const {
    return solid_.IsCellSolid(x, y);
}

bool Fluid2D::IsSolidWorld_(Real x, Real y) const {
    return solid_.IsWorldSolid(x, y);
}

bool Fluid2D::IsUFaceOpen_(int x, int y) const {
    return solid_.IsUFaceOpen(x, y);
}

bool Fluid2D::IsVFaceOpen_(int x, int y) const {
    return solid_.IsVFaceOpen(x, y);
}

Real Fluid2D::SolidU_(int x, int y) const {
    return solid_.UFaceVelocity(x, y);
}

Real Fluid2D::SolidV_(int x, int y) const {
    return solid_.VFaceVelocity(x, y);
}

Vector2R Fluid2D::ProjectOutOfSolid_(const Vector2R& position, const Vector2R& fallback) const {
    return solid_.ProjectOutOfSolid(position, fallback);
}

Real Fluid2D::SampleScalarWorld_(const CellCenteredScalarGrid2D& grid, Real x, Real y) const {
    // World-space cell center (i + 0.5) * dx maps to grid index i.
    const Real grid_x = CellCenterGridCoordinate2D(x, inverse_cell_size_);
    const Real grid_y = CellCenterGridCoordinate2D(y, inverse_cell_size_);
    return grid.Sample(grid_x, grid_y, config_.advection_interpolation);
}

Vector2R Fluid2D::SampleVelocityWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    return {SampleUWorld_(grid, x, y), SampleVWorld_(grid, x, y)};
}

Real Fluid2D::SampleUWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    // U is face-aligned in x and cell-centered in y.
    const Real grid_x = FaceGridCoordinate2D(x, inverse_cell_size_);
    const Real grid_y = CellCenterGridCoordinate2D(y, inverse_cell_size_);
    return grid.SampleU(grid_x, grid_y, config_.advection_interpolation);
}

Real Fluid2D::SampleVWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const {
    // V is cell-centered in x and face-aligned in y.
    const Real grid_x = CellCenterGridCoordinate2D(x, inverse_cell_size_);
    const Real grid_y = FaceGridCoordinate2D(y, inverse_cell_size_);
    return grid.SampleV(grid_x, grid_y, config_.advection_interpolation);
}

}  // namespace skyspaces
