#include "Fluid2D.h"
#include "Interpolation2D.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace skyspaces {

Fluid2D::Fluid2D(Config config) : config_(config) {
    InitializeFromConfig_();
}

void Fluid2D::Reset() {
    // density_.Fill(0.0);
    // density_tmp_.Fill(0.0);
    // pressure_.Fill(0.0);
    // divergence_.Fill(0.0);
    // vorticity_.Fill(0.0);
    velocity_.Fill(0.0);
    // velocity_x_tmp_.Fill(0.0);
    // velocity_y_tmp_.Fill(0.0);
    force_.Fill(0.0);
    // external_force_x_.Fill(0.0);
    // external_force_y_.Fill(0.0);
    // time_ = 0.0;
    // wind_phase_ = 0.0;
    // step_count_ = 0;
    // last_pressure_iterations_ = 0;
    // last_pressure_residual_ = 0.0;

    // if (_config.source_enabled) {
    //     ApplyConfiguredSource();
    // }
}

void Fluid2D::Step() {
    Step(config_.time_step);
}

void Fluid2D::Step(double dt) {
    if (dt <= 0.0) {
        return;
    }

    const int substeps = ComputeSubstepCount_(dt);
    const double substep_dt = dt / static_cast<double>(substeps);

    for (int substep = 0; substep < substeps; ++substep) {
        // if (config_.source_enabled) {
        //     ApplyConfiguredSource();
        // }

        force_.Fill(0.0);

        // ApplyBuoyancy();
        // if (_config.wind_enabled) {
        //     ApplyWind();
        // }

        // ComputeVorticityField();
        // ApplyVorticityConfinement();
        // ApplyAccumulatedForces(substep_dt);
        // DiffuseVelocity(substep_dt);
        // SolvePressureProjection(substep_dt);
        // AdvectFields(substep_dt);
        // ComputeVorticityField();

        time_ += substep_dt;
    }

    // ClearForces();
    ++step_count_;
}

void Fluid2D::InitializeFromConfig_() {
    assert(config_.resolution_x > 0 && config_.resolution_y > 0);

    inverse_cell_size_ = 1.0 / config_.cell_size;

    if (config_.time_step <= 0.0) {
        config_.time_step = std::max(1e-6, config_.cfl_number * config_.cell_size); // time_step <= cfl_number * cell_size / max_speed (max_speed=1.0)
    }

    // density_.Resize(_config.resolution_x, _config.resolution_y, 0.0);
    // density_tmp_.Resize(_config.resolution_x, _config.resolution_y, 0.0);
    // pressure_.Resize(_config.resolution_x, _config.resolution_y, 0.0);
    // divergence_.Resize(_config.resolution_x, _config.resolution_y, 0.0);
    // vorticity_.Resize(_config.resolution_x, _config.resolution_y, 0.0);
    velocity_.Resize(config_.resolution_x + 1, config_.resolution_y, 0.0);
    // velocity_x_tmp_.Resize(config_.resolution_x + 1, config_.resolution_y, 0.0);
    // velocity_y_tmp_.Resize(config_.resolution_x, config_.resolution_y + 1, 0.0);
    force_.Resize(config_.resolution_x + 1, config_.resolution_y, 0.0);
    // external_force_y_.Resize(_config.resolution_x, _config.resolution_y + 1, 0.0);
}

int Fluid2D::ComputeSubstepCount_(double dt) const {
    const double max_speed = ComputeMaxSpeed_();
    if (max_speed <= 0.0) {
        return 1;
    }

    // time_step <= cfl_number * cell_size / max_speed
    const double max_dt = config_.cfl_number * config_.cell_size / max_speed;
    const int count = static_cast<int>(std::ceil(dt / max_dt));
    return std::clamp(count, 1, config_.max_substeps);
}

double Fluid2D::ComputeMaxSpeed_() const {
    double max_speed = 0.0;
    for (double v : velocity_.UData()) {
        max_speed = std::max(max_speed, std::abs(v));
    }
    for (double v : velocity_.VData()) {
        max_speed = std::max(max_speed, std::abs(v));
    }
    return max_speed;
}

void Fluid2D::AdvectFields_(double dt) {
    SemiLagrangianBackwards_(dt, config_.semi_lagrangian_substeps);
    

    // for (int y = 0; y < config_.resolution_y; ++y) {
    //     for (int x = 0; x < config_.resolution_x; ++x) {
    //         const double nx = safeNormalizeCenter(x, config_.resolution_x);
    //         const double ny = safeNormalizeCenter(y, config_.resolution_y);
    //         const Vector2D velocity = SampleVelocityAtCenterCoords(nx, ny);
    //         const double back_x = clampValue(nx - velocity.x * dt / config_.domain_size_x, 0.0, 1.0);
    //         const double back_y = clampValue(ny - velocity.y * dt / config_.domain_size_y, 0.0, 1.0);
    //         density_tmp_(x, y) = SampleDensityAtCenterCoords(back_x, back_y);
    //     }
    // }
    // density_.Data() = density_tmp_.Data();

    // for (int y = 0; y < velocity_x_.Height(); ++y) {
    //     for (int x = 0; x < velocity_x_.Width(); ++x) {
    //         const double nx = safeNormalizeX(x, _config.resolution_x);
    //         const double ny = safeNormalizeCenter(y, _config.resolution_y);
    //         const Vector2D velocity = SampleVelocityAtUFaceCoords(nx, ny);
    //         const double back_x = clampValue(nx - velocity.x * dt / _config.domain_size_x, 0.0, 1.0);
    //         const double back_y = clampValue(ny - velocity.y * dt / _config.domain_size_y, 0.0, 1.0);
    //         velocity_x_tmp_(x, y) = SampleUAtCenterCoords(back_x, back_y);
    //     }
    // }
    // velocity_x_.Data() = velocity_x_tmp_.Data();

    // for (int y = 0; y < velocity_y_.Height(); ++y) {
    //     for (int x = 0; x < velocity_y_.Width(); ++x) {
    //         const double nx = safeNormalizeCenter(x, _config.resolution_x);
    //         const double ny = safeNormalizeX(y, _config.resolution_y);
    //         const Vector2D velocity = SampleVelocityAtVFaceCoords(nx, ny);
    //         const double back_x = clampValue(nx - velocity.x * dt / _config.domain_size_x, 0.0, 1.0);
    //         const double back_y = clampValue(ny - velocity.y * dt / _config.domain_size_y, 0.0, 1.0);
    //         velocity_y_tmp_(x, y) = SampleVAtCenterCoords(back_x, back_y);
    //     }
    // }
    // velocity_y_.Data() = velocity_y_tmp_.Data();

    // SetBoundaryConditions();
}

void Fluid2D::SemiLagrangianBackwards_(double dt, int substeps) {
    const double substep_dt = dt / static_cast<double>(substeps);
    
    int width = config_.resolution_x;
    int height = config_.resolution_y;
    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            semi_lagrangian_position_.U(i, j) = static_cast<double>(i) + 0.5;
            semi_lagrangian_position_.V(i, j) = static_cast<double>(j) + 0.5;
        }
    }
    
    for (int step = 0; step < substeps; ++step) {
        for (int j = 0; j < height; ++j) {
            for (int i = 0; i < width; ++i) {
                const double x = semi_lagrangian_position_.U(i, j);
                const double y = semi_lagrangian_position_.V(i, j);
                const Vector2D velocity = velocity_.Sample(x, y);
                semi_lagrangian_position_.U(i, j) = x - velocity.x * substep_dt * inverse_cell_size_;
                semi_lagrangian_position_.V(i, j) = y - velocity.y * substep_dt * inverse_cell_size_;
            }
        }
    }
}

// namespace {

// double clampValue(double value, double min_value, double max_value) {
//     return std::max(min_value, std::min(value, max_value));
// }

// int clampIndex(int value, int min_value, int max_value) {
//     return std::max(min_value, std::min(value, max_value));
// }

// std::pair<int, int> normalizedIndexRange(
//     double min_normalized,
//     double max_normalized,
//     int resolution) {
//     const double start_normalized = clampValue(std::min(min_normalized, max_normalized), 0.0, 1.0);
//     const double end_normalized = clampValue(std::max(min_normalized, max_normalized), 0.0, 1.0);

//     const int start_index = clampIndex(
//         static_cast<int>(std::floor(start_normalized * resolution)), 0, resolution - 1);
//     const int end_index = clampIndex(
//         static_cast<int>(std::ceil(end_normalized * resolution)), 0, resolution);

//     return {start_index, end_index};
// }

// double safeNormalizeX(int x, int resolution_x) {
//     return static_cast<double>(x) / static_cast<double>(std::max(resolution_x, 1));
// }

// double safeNormalizeCenter(int index, int resolution) {
//     return (static_cast<double>(index) + 0.5) / static_cast<double>(std::max(resolution, 1));
// }

// double sampleClamped(const MACScalarGrid2D& grid, double x, double y) {
//     return grid.Sample(clampValue(x, 0.0, static_cast<double>(grid.Width() - 1)),
//                        clampValue(y, 0.0, static_cast<double>(grid.Height() - 1)));
// }

// }  // namespace


// void Fluid2D::AddDensitySource(
//     double min_x,
//     double max_x,
//     double min_y,
//     double max_y,
//     double density) {
//     const auto [start_x, end_x] = normalizedIndexRange(min_x, max_x, _config.resolution_x);
//     const auto [start_y, end_y] = normalizedIndexRange(min_y, max_y, _config.resolution_y);

//     for (int y = start_y; y < end_y; ++y) {
//         for (int x = start_x; x < end_x; ++x) {
//             density_(x, y) = std::max(density_(x, y), density);
//         }
//     }
// }

// void Fluid2D::AddForce(
//     double min_x,
//     double max_x,
//     double min_y,
//     double max_y,
//     double force_x,
//     double force_y) {
//     const double x0 = clampValue(std::min(min_x, max_x), 0.0, 1.0);
//     const double x1 = clampValue(std::max(min_x, max_x), 0.0, 1.0);
//     const double y0 = clampValue(std::min(min_y, max_y), 0.0, 1.0);
//     const double y1 = clampValue(std::max(min_y, max_y), 0.0, 1.0);

//     for (int y = 0; y < external_force_x_.Height(); ++y) {
//         const double py = safeNormalizeCenter(y, _config.resolution_y);
//         if (py < y0 || py > y1) {
//             continue;
//         }

//         for (int x = 0; x < external_force_x_.Width(); ++x) {
//             const double px = safeNormalizeX(x, _config.resolution_x);
//             if (px >= x0 && px <= x1) {
//                 external_force_x_(x, y) += force_x;
//             }
//         }
//     }

//     for (int y = 0; y < external_force_y_.Height(); ++y) {
//         const double py = safeNormalizeX(y, _config.resolution_y);
//         if (py < y0 || py > y1) {
//             continue;
//         }

//         for (int x = 0; x < external_force_y_.Width(); ++x) {
//             const double px = safeNormalizeCenter(x, _config.resolution_x);
//             if (px >= x0 && px <= x1) {
//                 external_force_y_(x, y) += force_y;
//             }
//         }
//     }
// }

// void Fluid2D::ClearDensity() {
//     density_.Fill(0.0);
//     density_tmp_.Fill(0.0);
// }

// void Fluid2D::ClearForces() {
//     external_force_x_.Fill(0.0);
//     external_force_y_.Fill(0.0);
//     force_x_.Fill(0.0);
//     force_y_.Fill(0.0);
// }

// int Fluid2D::ResolutionX() const noexcept {
//     return _config.resolution_x;
// }

// int Fluid2D::ResolutionY() const noexcept {
//     return _config.resolution_y;
// }

// double Fluid2D::CellSize() const noexcept {
//     return cell_size_;
// }

// double Fluid2D::InverseCellSize() const noexcept {
//     return inverse_cell_size_;
// }

// double Fluid2D::DomainSizeX() const noexcept {
//     return _config.domain_size_x;
// }

// double Fluid2D::DomainSizeY() const noexcept {
//     return _config.domain_size_y;
// }

// void Fluid2D::ApplyConfiguredSource() {
//     AddDensitySource(
//         _config.source_min_x,
//         _config.source_max_x,
//         _config.source_min_y,
//         _config.source_max_y,
//         _config.source_density);

//     AddForce(
//         _config.source_min_x,
//         _config.source_max_x,
//         _config.source_min_y,
//         _config.source_max_y,
//         0.0,
//         _config.buoyancy);
// }

// void Fluid2D::ApplyBuoyancy() {
//     for (int y = 0; y < velocity_y_.Height(); ++y) {
//         const double py = safeNormalizeX(y, _config.resolution_y);
//         for (int x = 0; x < velocity_y_.Width(); ++x) {
//             const double px = safeNormalizeCenter(x, _config.resolution_x);
//             const double density_sample = SampleDensityAtCenterCoords(px, py);
//             force_y_(x, y) += _config.buoyancy * density_sample;
//         }
//     }
// }

// void Fluid2D::ApplyWind() {
//     for (int y = 0; y < velocity_x_.Height(); ++y) {
//         for (int x = 0; x < velocity_x_.Width(); ++x) {
//             force_x_(x, y) += _config.wind_strength;
//         }
//     }
// }

// void Fluid2D::ApplyVorticityConfinement() {
//     const int w = _config.resolution_x;
//     const int h = _config.resolution_y;

//     for (int y = 1; y < h - 1; ++y) {
//         for (int x = 1; x < w - 1; ++x) {
//             const double dwdx = (std::abs(vorticity_(x + 1, y)) - std::abs(vorticity_(x - 1, y))) * 0.5 * inverse_cell_size_;
//             const double dwdy = (std::abs(vorticity_(x, y + 1)) - std::abs(vorticity_(x, y - 1))) * 0.5 * inverse_cell_size_;
//             const double length = std::sqrt(dwdx * dwdx + dwdy * dwdy) + 1e-9;
//             const double nx = dwdx / length;
//             const double ny = dwdy / length;
//             const double vortex = vorticity_(x, y);

//             const double force_x = _config.vorticity_confinement * (ny * vortex);
//             const double force_y = _config.vorticity_confinement * (-nx * vortex);

//             const int iu = std::clamp(x, 0, velocity_x_.Width() - 1);
//             const int iv = std::clamp(x, 0, velocity_y_.Width() - 1);
//             const int ju = std::clamp(y, 0, velocity_x_.Height() - 1);
//             const int jv = std::clamp(y, 0, velocity_y_.Height() - 1);

//             force_x_(iu, ju) += force_x;
//             force_y_(iv, jv) += force_y;
//         }
//     }
// }

// void Fluid2D::ApplyAccumulatedForces(double dt) {
//     const double scale = dt;

//     for (int y = 0; y < velocity_x_.Height(); ++y) {
//         for (int x = 0; x < velocity_x_.Width(); ++x) {
//             velocity_x_(x, y) += scale * (force_x_(x, y) + external_force_x_(x, y));
//         }
//     }

//     for (int y = 0; y < velocity_y_.Height(); ++y) {
//         for (int x = 0; x < velocity_y_.Width(); ++x) {
//             velocity_y_(x, y) += scale * (force_y_(x, y) + external_force_y_(x, y));
//         }
//     }
// }

// void Fluid2D::DiffuseVelocity(double dt) {
//     if (_config.viscosity <= 0.0) {
//         return;
//     }

//     const double alpha = _config.viscosity * dt * inverse_cell_size_ * inverse_cell_size_;
//     const double inv_coefficient = 1.0 / (1.0 + 4.0 * alpha);

//     auto diffuseGrid = [&](MACScalarGrid2D& field, MACScalarGrid2D& scratch) {
//         for (int iter = 0; iter < _config.diffusion_iterations; ++iter) {
//             for (int y = 1; y < field.Height() - 1; ++y) {
//                 for (int x = 1; x < field.Width() - 1; ++x) {
//                     const double sum_neighbors = field(x - 1, y) + field(x + 1, y) + field(x, y - 1) + field(x, y + 1);
//                     scratch(x, y) = (field(x, y) + alpha * sum_neighbors) * inv_coefficient;
//                 }
//             }

//             for (int y = 1; y < field.Height() - 1; ++y) {
//                 for (int x = 1; x < field.Width() - 1; ++x) {
//                     field(x, y) = scratch(x, y);
//                 }
//             }
//             SetBoundaryConditions();
//         }
//     };

//     diffuseGrid(velocity_x_, velocity_x_tmp_);
//     diffuseGrid(velocity_y_, velocity_y_tmp_);
// }

// void Fluid2D::SolvePressureProjection(double dt) {
//     SetBoundaryConditions();
//     ComputeDivergenceField();

//     pressure_.Fill(0.0);
//     SolvePoisson(dt);
//     CorrectVelocity(dt);
//     SetBoundaryConditions();
// }

// void Fluid2D::SetBoundaryConditions() {
//     const int ux = velocity_x_.Width();
//     const int uy = velocity_x_.Height();
//     for (int y = 0; y < uy; ++y) {
//         velocity_x_(0, y) = 0.0;
//         velocity_x_(ux - 1, y) = 0.0;
//     }

//     const int vx = velocity_y_.Width();
//     const int vy = velocity_y_.Height();
//     for (int x = 0; x < vx; ++x) {
//         velocity_y_(x, 0) = 0.0;
//         velocity_y_(x, vy - 1) = 0.0;
//     }

//     for (int x = 0; x < pressure_.Width(); ++x) {
//         pressure_(x, 0) = 0.0;
//         pressure_(x, pressure_.Height() - 1) = 0.0;
//     }
//     for (int y = 0; y < pressure_.Height(); ++y) {
//         pressure_(0, y) = 0.0;
//         pressure_(pressure_.Width() - 1, y) = 0.0;
//     }
// }

// void Fluid2D::ComputeDivergenceField() {
//     for (int y = 0; y < _config.resolution_y; ++y) {
//         for (int x = 0; x < _config.resolution_x; ++x) {
//             const double du = velocity_x_(x + 1, y) - velocity_x_(x, y);
//             const double dv = velocity_y_(x, y + 1) - velocity_y_(x, y);
//             divergence_(x, y) = (du + dv) * inverse_cell_size_;
//         }
//     }
// }

// void Fluid2D::SolvePoisson(double dt) {
//     const double alpha = cell_size_ * cell_size_ / dt;
//     const int w = pressure_.Width();
//     const int h = pressure_.Height();

//     for (int iter = 0; iter < _config.pressure_iterations; ++iter) {
//         double max_change = 0.0;

//         for (int y = 1; y < h - 1; ++y) {
//             for (int x = 1; x < w - 1; ++x) {
//                 const double sum_neighbors = pressure_(x - 1, y) + pressure_(x + 1, y) + pressure_(x, y - 1) + pressure_(x, y + 1);
//                 const double new_pressure = (sum_neighbors - alpha * divergence_(x, y)) * 0.25;
//                 const double change = std::abs(new_pressure - pressure_(x, y));
//                 pressure_(x, y) = new_pressure;
//                 max_change = std::max(max_change, change);
//             }
//         }

//         SetBoundaryConditions();
//         if (max_change < _config.pressure_tolerance) {
//             last_pressure_iterations_ = iter + 1;
//             last_pressure_residual_ = max_change;
//             return;
//         }
//     }

//     last_pressure_iterations_ = _config.pressure_iterations;
//     last_pressure_residual_ = 0.0;
// }

// void Fluid2D::CorrectVelocity(double dt) {
//     const double scale = dt * inverse_cell_size_;

//     for (int y = 0; y < velocity_x_.Height(); ++y) {
//         for (int x = 1; x < velocity_x_.Width() - 1; ++x) {
//             const double p_left = PressureAtWithNeumann(x - 1, y);
//             const double p_right = PressureAtWithNeumann(x, y);
//             velocity_x_(x, y) -= scale * (p_right - p_left);
//         }
//     }

//     for (int y = 1; y < velocity_y_.Height() - 1; ++y) {
//         for (int x = 0; x < velocity_y_.Width(); ++x) {
//             const double p_down = PressureAtWithNeumann(x, y - 1);
//             const double p_up = PressureAtWithNeumann(x, y);
//             velocity_y_(x, y) -= scale * (p_up - p_down);
//         }
//     }
// }

// void Fluid2D::ComputeVorticityField() {
//     for (int y = 0; y < _config.resolution_y; ++y) {
//         for (int x = 0; x < _config.resolution_x; ++x) {
//             if (x == _config.resolution_x - 1 || y == _config.resolution_y - 1) {
//                 vorticity_(x, y) = 0.0;
//                 continue;
//             }

//             const double du_dy = (velocity_x_(x, y + 1) - velocity_x_(x, y)) * inverse_cell_size_;
//             const double dv_dx = (velocity_y_(x + 1, y) - velocity_y_(x, y)) * inverse_cell_size_;
//             vorticity_(x, y) = dv_dx - du_dy;
//         }
//     }
// }







// double Fluid2D::SampleDensityAtCenterCoords(double x, double y) const {
//     const double grid_x = x * static_cast<double>(_config.resolution_x) - 0.5;
//     const double grid_y = y * static_cast<double>(_config.resolution_y) - 0.5;
//     return sampleClamped(density_, grid_x, grid_y);
// }

// double Fluid2D::SampleUAtCenterCoords(double x, double y) const {
//     const double grid_x = x * static_cast<double>(_config.resolution_x);
//     const double grid_y = y * static_cast<double>(_config.resolution_y) - 0.5;
//     return sampleClamped(velocity_x_, grid_x, grid_y);
// }

// double Fluid2D::SampleVAtCenterCoords(double x, double y) const {
//     const double grid_x = x * static_cast<double>(_config.resolution_x) - 0.5;
//     const double grid_y = y * static_cast<double>(_config.resolution_y);
//     return sampleClamped(velocity_y_, grid_x, grid_y);
// }

// Vector2D Fluid2D::SampleVelocityAtCenterCoords(double x, double y) const {
//     return {
//         SampleUAtCenterCoords(x, y),
//         SampleVAtCenterCoords(x, y)
//     };
// }

// Vector2D Fluid2D::SampleVelocityAtUFaceCoords(double x, double y) const {
//     return {SampleUAtCenterCoords(x, y), SampleVAtCenterCoords(x, y)};
// }

// Vector2D Fluid2D::SampleVelocityAtVFaceCoords(double x, double y) const {
//     return {SampleUAtCenterCoords(x, y), SampleVAtCenterCoords(x, y)};
// }

// double Fluid2D::PressureAtWithNeumann(int x, int y) const {
//     const int clamped_x = std::clamp(x, 0, pressure_.Width() - 1);
//     const int clamped_y = std::clamp(y, 0, pressure_.Height() - 1);
//     return pressure_(clamped_x, clamped_y);
// }

// double Fluid2D::TimeStep() const noexcept {
//     return _config.time_step;
// }

// double Fluid2D::Time() const noexcept {
//     return time_;
// }

// int Fluid2D::StepCount() const noexcept {
//     return step_count_;
// }

// int Fluid2D::LastPressureIterations() const noexcept {
//     return last_pressure_iterations_;
// }

// double Fluid2D::LastPressureResidual() const noexcept {
//     return last_pressure_residual_;
// }

// Vector2D Fluid2D::SampleVelocityNormalized(double x, double y) const {
//     const double clamped_x = clampValue(x, 0.0, 1.0);
//     const double clamped_y = clampValue(y, 0.0, 1.0);
//     return SampleVelocityAtCenterCoords(clamped_x, clamped_y);
// }

// double Fluid2D::SampleDensityNormalized(double x, double y) const {
//     const double clamped_x = clampValue(x, 0.0, 1.0);
//     const double clamped_y = clampValue(y, 0.0, 1.0);
//     return SampleDensityAtCenterCoords(clamped_x, clamped_y);
// }

// const Fluid2D::Config& Fluid2D::Configs() const noexcept {
//     return _config;
// }

// Fluid2D::Config& Fluid2D::Configs() noexcept {
//     return _config;
// }

// const MACScalarGrid2D& Fluid2D::Density() const noexcept {
//     return density_;
// }

// const MACScalarGrid2D& Fluid2D::Pressure() const noexcept {
//     return pressure_;
// }

// const MACScalarGrid2D& Fluid2D::Divergence() const noexcept {
//     return divergence_;
// }

// const MACScalarGrid2D& Fluid2D::Vorticity() const noexcept {
//     return vorticity_;
// }

// const MACScalarGrid2D& Fluid2D::VelocityX() const noexcept {
//     return velocity_x_;
// }

// const MACScalarGrid2D& Fluid2D::VelocityY() const noexcept {
//     return velocity_y_;
// }

}  // namespace skyspaces
