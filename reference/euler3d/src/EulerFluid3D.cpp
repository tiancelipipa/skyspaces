#include "EulerFluid3D.h"

#include <Eigen/SparseCore>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace skyspaces::reference {
namespace {

using Clock = std::chrono::steady_clock;

Real SecondsSince(Clock::time_point start) {
    return std::chrono::duration<Real>(Clock::now() - start).count();
}

int FlattenCellIndex3D(int x, int y, int z, int height, int depth) {
    return (x * height + y) * depth + z;
}

Real CellCenterGridCoordinate(Real world, Real inverse_cell_size) {
    return world * inverse_cell_size - 0.5;
}

Real FaceGridCoordinate(Real world, Real inverse_cell_size) {
    return world * inverse_cell_size;
}

Real MaxAbs(const std::vector<Real>& values) {
    Real result = 0.0;
    for (const Real value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

}  // namespace

CellCenteredScalarGrid3D::CellCenteredScalarGrid3D(int width, int height, int depth, Real value) {
    Resize(width, height, depth, value);
}

void CellCenteredScalarGrid3D::Resize(int width, int height, int depth, Real value) {
    assert(width >= 0 && height >= 0 && depth >= 0);
    width_ = width;
    height_ = height;
    depth_ = depth;
    data_.assign(
        static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_) * static_cast<std::size_t>(depth_),
        value);
}

void CellCenteredScalarGrid3D::Fill(Real value) {
    std::fill(data_.begin(), data_.end(), value);
}

int CellCenteredScalarGrid3D::Width() const noexcept {
    return width_;
}

int CellCenteredScalarGrid3D::Height() const noexcept {
    return height_;
}

int CellCenteredScalarGrid3D::Depth() const noexcept {
    return depth_;
}

bool CellCenteredScalarGrid3D::Empty() const noexcept {
    return data_.empty();
}

Real& CellCenteredScalarGrid3D::operator()(int x, int y, int z) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z < depth_);
    return data_[static_cast<std::size_t>(Index_(x, y, z))];
}

Real CellCenteredScalarGrid3D::operator()(int x, int y, int z) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z < depth_);
    return data_[static_cast<std::size_t>(Index_(x, y, z))];
}

Real CellCenteredScalarGrid3D::Sample(Real x, Real y, Real z) const {
    return Sample(x, y, z, InterpolationMethod3D::Trilinear);
}

Real CellCenteredScalarGrid3D::Sample(Real x, Real y, Real z, InterpolationMethod3D method) const {
    return skyspaces::reference::Sample(data_, width_, height_, depth_, x, y, z, method);
}

const std::vector<Real>& CellCenteredScalarGrid3D::Data() const noexcept {
    return data_;
}

std::vector<Real>& CellCenteredScalarGrid3D::Data() noexcept {
    return data_;
}

int CellCenteredScalarGrid3D::Index_(int x, int y, int z) const noexcept {
    return FlattenCellIndex3D(x, y, z, height_, depth_);
}

FaceCenteredVectorGrid3D::FaceCenteredVectorGrid3D(int width, int height, int depth, Real value) {
    Resize(width, height, depth, value);
}

void FaceCenteredVectorGrid3D::Resize(int width, int height, int depth, Real value) {
    assert(width >= 0 && height >= 0 && depth >= 0);
    width_ = width;
    height_ = height;
    depth_ = depth;
    u_data_.assign(
        static_cast<std::size_t>(width_ + 1) * static_cast<std::size_t>(height_) *
            static_cast<std::size_t>(depth_),
        value);
    v_data_.assign(
        static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_ + 1) *
            static_cast<std::size_t>(depth_),
        value);
    w_data_.assign(
        static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_) *
            static_cast<std::size_t>(depth_ + 1),
        value);
}

void FaceCenteredVectorGrid3D::Fill(Real value) {
    std::fill(u_data_.begin(), u_data_.end(), value);
    std::fill(v_data_.begin(), v_data_.end(), value);
    std::fill(w_data_.begin(), w_data_.end(), value);
}

int FaceCenteredVectorGrid3D::Width() const noexcept {
    return width_;
}

int FaceCenteredVectorGrid3D::Height() const noexcept {
    return height_;
}

int FaceCenteredVectorGrid3D::Depth() const noexcept {
    return depth_;
}

int FaceCenteredVectorGrid3D::UWidth() const noexcept {
    return width_ + 1;
}

int FaceCenteredVectorGrid3D::UHeight() const noexcept {
    return height_;
}

int FaceCenteredVectorGrid3D::UDepth() const noexcept {
    return depth_;
}

int FaceCenteredVectorGrid3D::VWidth() const noexcept {
    return width_;
}

int FaceCenteredVectorGrid3D::VHeight() const noexcept {
    return height_ + 1;
}

int FaceCenteredVectorGrid3D::VDepth() const noexcept {
    return depth_;
}

int FaceCenteredVectorGrid3D::WWidth() const noexcept {
    return width_;
}

int FaceCenteredVectorGrid3D::WHeight() const noexcept {
    return height_;
}

int FaceCenteredVectorGrid3D::WDepth() const noexcept {
    return depth_ + 1;
}

Real& FaceCenteredVectorGrid3D::U(int x, int y, int z) {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z < depth_);
    return u_data_[static_cast<std::size_t>(UIndex_(x, y, z))];
}

Real FaceCenteredVectorGrid3D::U(int x, int y, int z) const {
    assert(x >= 0 && x <= width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z < depth_);
    return u_data_[static_cast<std::size_t>(UIndex_(x, y, z))];
}

Real& FaceCenteredVectorGrid3D::V(int x, int y, int z) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    assert(z >= 0 && z < depth_);
    return v_data_[static_cast<std::size_t>(VIndex_(x, y, z))];
}

Real FaceCenteredVectorGrid3D::V(int x, int y, int z) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y <= height_);
    assert(z >= 0 && z < depth_);
    return v_data_[static_cast<std::size_t>(VIndex_(x, y, z))];
}

Real& FaceCenteredVectorGrid3D::W(int x, int y, int z) {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z <= depth_);
    return w_data_[static_cast<std::size_t>(WIndex_(x, y, z))];
}

Real FaceCenteredVectorGrid3D::W(int x, int y, int z) const {
    assert(x >= 0 && x < width_);
    assert(y >= 0 && y < height_);
    assert(z >= 0 && z <= depth_);
    return w_data_[static_cast<std::size_t>(WIndex_(x, y, z))];
}

Real FaceCenteredVectorGrid3D::SampleU(Real x, Real y, Real z) const {
    return SampleU(x, y, z, InterpolationMethod3D::Trilinear);
}

Real FaceCenteredVectorGrid3D::SampleU(Real x, Real y, Real z, InterpolationMethod3D method) const {
    return skyspaces::reference::Sample(u_data_, UWidth(), UHeight(), UDepth(), x, y, z, method);
}

Real FaceCenteredVectorGrid3D::SampleV(Real x, Real y, Real z) const {
    return SampleV(x, y, z, InterpolationMethod3D::Trilinear);
}

Real FaceCenteredVectorGrid3D::SampleV(Real x, Real y, Real z, InterpolationMethod3D method) const {
    return skyspaces::reference::Sample(v_data_, VWidth(), VHeight(), VDepth(), x, y, z, method);
}

Real FaceCenteredVectorGrid3D::SampleW(Real x, Real y, Real z) const {
    return SampleW(x, y, z, InterpolationMethod3D::Trilinear);
}

Real FaceCenteredVectorGrid3D::SampleW(Real x, Real y, Real z, InterpolationMethod3D method) const {
    return skyspaces::reference::Sample(w_data_, WWidth(), WHeight(), WDepth(), x, y, z, method);
}

Vector3R FaceCenteredVectorGrid3D::Sample(Real x, Real y, Real z) const {
    return Sample(x, y, z, InterpolationMethod3D::Trilinear);
}

Vector3R FaceCenteredVectorGrid3D::Sample(Real x, Real y, Real z, InterpolationMethod3D method) const {
    return {SampleU(x, y, z, method), SampleV(x, y, z, method), SampleW(x, y, z, method)};
}

const std::vector<Real>& FaceCenteredVectorGrid3D::UData() const noexcept {
    return u_data_;
}

std::vector<Real>& FaceCenteredVectorGrid3D::UData() noexcept {
    return u_data_;
}

const std::vector<Real>& FaceCenteredVectorGrid3D::VData() const noexcept {
    return v_data_;
}

std::vector<Real>& FaceCenteredVectorGrid3D::VData() noexcept {
    return v_data_;
}

const std::vector<Real>& FaceCenteredVectorGrid3D::WData() const noexcept {
    return w_data_;
}

std::vector<Real>& FaceCenteredVectorGrid3D::WData() noexcept {
    return w_data_;
}

int FaceCenteredVectorGrid3D::UIndex_(int x, int y, int z) const noexcept {
    return FlattenCellIndex3D(x, y, z, height_, depth_);
}

int FaceCenteredVectorGrid3D::VIndex_(int x, int y, int z) const noexcept {
    return (x * (height_ + 1) + y) * depth_ + z;
}

int FaceCenteredVectorGrid3D::WIndex_(int x, int y, int z) const noexcept {
    return (x * height_ + y) * (depth_ + 1) + z;
}

EulerFluid3D::EulerFluid3D(Config config) : config_(config) {
    InitializeFromConfig_();
}

void EulerFluid3D::Reset() {
    smoke_density_.Fill(0.0);
    smoke_density_tmp_.Fill(0.0);
    temperature_.Fill(config_.ambient_temperature);
    temperature_tmp_.Fill(config_.ambient_temperature);
    pressure_.Fill(0.0);
    divergence_.Fill(0.0);
    velocity_.Fill(0.0);
    velocity_tmp_.Fill(0.0);
    ConfigureAdvectionScratch_();
    time_ = 0.0;
    step_count_ = 0;
    last_pressure_iterations_ = 0;
    last_pressure_residual_ = 0.0;
    last_step_timings_ = {};
}

void EulerFluid3D::Step() {
    Step(config_.time_step);
}

void EulerFluid3D::Step(Real dt) {
    if (dt <= 0.0) {
        return;
    }

    last_step_timings_ = {};
    const auto step_start = Clock::now();
    const int substeps = ComputeSubstepCount_(dt);
    last_step_timings_.substep_count = substeps;
    const Real substep_dt = dt / static_cast<Real>(substeps);
    for (int substep = 0; substep < substeps; ++substep) {
        auto stage_start = Clock::now();
        AdvectVelocity_(substep_dt);
        last_step_timings_.advect_velocity += SecondsSince(stage_start);

        stage_start = Clock::now();
        ApplyConfiguredVelocitySource_(substep_dt);
        last_step_timings_.velocity_source += SecondsSince(stage_start);

        stage_start = Clock::now();
        SetBoundaryConditions_();
        const Real pre_projection_boundary_seconds = SecondsSince(stage_start);
        last_step_timings_.pre_projection_boundary += pre_projection_boundary_seconds;
        last_step_timings_.boundary_conditions += pre_projection_boundary_seconds;

        stage_start = Clock::now();
        ProjectVelocity_(substep_dt);
        last_step_timings_.projection += SecondsSince(stage_start);

        stage_start = Clock::now();
        AdvectScalars_(substep_dt);
        last_step_timings_.advect_scalars += SecondsSince(stage_start);

        stage_start = Clock::now();
        ApplyConfiguredScalarSource_(substep_dt);
        last_step_timings_.scalar_source += SecondsSince(stage_start);

        time_ += substep_dt;
    }

    last_step_timings_.total = SecondsSince(step_start);
    ++step_count_;
}

int EulerFluid3D::ResolutionX() const noexcept {
    return config_.resolution_x;
}

int EulerFluid3D::ResolutionY() const noexcept {
    return config_.resolution_y;
}

int EulerFluid3D::ResolutionZ() const noexcept {
    return config_.resolution_z;
}

Real EulerFluid3D::CellSize() const noexcept {
    return config_.cell_size;
}

Real EulerFluid3D::TimeStep() const noexcept {
    return config_.time_step;
}

Real EulerFluid3D::Time() const noexcept {
    return time_;
}

int EulerFluid3D::StepCount() const noexcept {
    return step_count_;
}

int EulerFluid3D::LastPressureIterations() const noexcept {
    return last_pressure_iterations_;
}

Real EulerFluid3D::LastPressureResidual() const noexcept {
    return last_pressure_residual_;
}

const EulerFluid3D::StepTimings& EulerFluid3D::LastStepTimings() const noexcept {
    return last_step_timings_;
}

Real EulerFluid3D::SampleSmokeNormalized(Real x, Real y, Real z) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    const Real world_z = std::clamp(z, 0.0, 1.0) * static_cast<Real>(config_.resolution_z) * config_.cell_size;
    return SampleScalarWorld_(smoke_density_, world_x, world_y, world_z);
}

Real EulerFluid3D::SampleTemperatureNormalized(Real x, Real y, Real z) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    const Real world_z = std::clamp(z, 0.0, 1.0) * static_cast<Real>(config_.resolution_z) * config_.cell_size;
    return SampleScalarWorld_(temperature_, world_x, world_y, world_z);
}

Vector3R EulerFluid3D::SampleVelocityNormalized(Real x, Real y, Real z) const {
    const Real world_x = std::clamp(x, 0.0, 1.0) * static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real world_y = std::clamp(y, 0.0, 1.0) * static_cast<Real>(config_.resolution_y) * config_.cell_size;
    const Real world_z = std::clamp(z, 0.0, 1.0) * static_cast<Real>(config_.resolution_z) * config_.cell_size;
    return SampleVelocityWorld_(velocity_, world_x, world_y, world_z);
}

const EulerFluid3D::Config& EulerFluid3D::Configs() const noexcept {
    return config_;
}

EulerFluid3D::Config& EulerFluid3D::Configs() noexcept {
    return config_;
}

const CellCenteredScalarGrid3D& EulerFluid3D::SmokeDensity() const noexcept {
    return smoke_density_;
}

CellCenteredScalarGrid3D& EulerFluid3D::SmokeDensity() noexcept {
    return smoke_density_;
}

const CellCenteredScalarGrid3D& EulerFluid3D::Temperature() const noexcept {
    return temperature_;
}

CellCenteredScalarGrid3D& EulerFluid3D::Temperature() noexcept {
    return temperature_;
}

const CellCenteredScalarGrid3D& EulerFluid3D::Pressure() const noexcept {
    return pressure_;
}

const CellCenteredScalarGrid3D& EulerFluid3D::Divergence() const noexcept {
    return divergence_;
}

const FaceCenteredVectorGrid3D& EulerFluid3D::Velocity() const noexcept {
    return velocity_;
}

void EulerFluid3D::InitializeFromConfig_() {
    assert(config_.numeric_epsilon > 0.0);
    assert(config_.resolution_x > 0 && config_.resolution_y > 0 && config_.resolution_z > 0);
    assert(config_.cell_size > 0.0);
    assert(config_.time_step > 0.0);
    assert(config_.max_substeps > 0);
    assert(config_.cfl_number >= 1.0 && config_.cfl_number <= 10.0);
    assert(config_.fluid_density > 0.0);
    assert(config_.pressure_iterations > 0);
    assert(config_.pressure_tolerance > 0.0);
    assert(config_.max_smoke_density > 0.0);
    assert(config_.source_min_x >= 0.0 && config_.source_min_x < config_.source_max_x && config_.source_max_x <= 1.0);
    assert(config_.source_min_y >= 0.0 && config_.source_min_y < config_.source_max_y && config_.source_max_y <= 1.0);
    assert(config_.source_min_z >= 0.0 && config_.source_min_z < config_.source_max_z && config_.source_max_z <= 1.0);
    assert(config_.source_temperature_rate >= 0.0);

    inverse_cell_size_ = 1.0 / config_.cell_size;
    smoke_density_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    smoke_density_tmp_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    temperature_.Resize(
        config_.resolution_x,
        config_.resolution_y,
        config_.resolution_z,
        config_.ambient_temperature);
    temperature_tmp_.Resize(
        config_.resolution_x,
        config_.resolution_y,
        config_.resolution_z,
        config_.ambient_temperature);
    pressure_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    divergence_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    velocity_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    velocity_tmp_.Resize(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    ConfigureAdvectionScratch_();
}

void EulerFluid3D::ConfigureAdvectionScratch_() {
    if (config_.advection_scheme == AdvectionScheme3D::SemiLagrangian) {
        ClearAdvectionScratch_();
        return;
    }

    velocity_first_pass_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    velocity_back_pass_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    velocity_corrected_source_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    smoke_first_pass_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    smoke_back_pass_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    smoke_corrected_source_.emplace(config_.resolution_x, config_.resolution_y, config_.resolution_z, 0.0);
    temperature_first_pass_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.resolution_z,
        config_.ambient_temperature);
    temperature_back_pass_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.resolution_z,
        config_.ambient_temperature);
    temperature_corrected_source_.emplace(
        config_.resolution_x,
        config_.resolution_y,
        config_.resolution_z,
        config_.ambient_temperature);
}

void EulerFluid3D::ClearAdvectionScratch_() {
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

int EulerFluid3D::ComputeSubstepCount_(Real dt) const {
    const Real max_speed = ComputeMaxSpeed_();
    if (max_speed <= config_.numeric_epsilon) {
        return 1;
    }

    const Real max_dt = config_.cfl_number * config_.cell_size / max_speed;
    const int count = static_cast<int>(std::ceil(dt / std::max(max_dt, config_.numeric_epsilon)));
    return std::clamp(count, 1, config_.max_substeps);
}

Real EulerFluid3D::ComputeMaxSpeed_() const {
    return std::max({MaxAbs(velocity_.UData()), MaxAbs(velocity_.VData()), MaxAbs(velocity_.WData())});
}

void EulerFluid3D::AdvectVelocity_(Real dt) {
    switch (config_.advection_scheme) {
        case AdvectionScheme3D::SemiLagrangian:
            AdvectVelocitySemiLagrangian3D(*this, dt);
            break;
        case AdvectionScheme3D::MacCormackBFECC:
            AdvectVelocityMacCormackBFECC3D(*this, dt);
            break;
    }

    velocity_.UData() = velocity_tmp_.UData();
    velocity_.VData() = velocity_tmp_.VData();
    velocity_.WData() = velocity_tmp_.WData();
}

void EulerFluid3D::AdvectScalars_(Real dt) {
    auto stage_start = Clock::now();
    switch (config_.advection_scheme) {
        case AdvectionScheme3D::SemiLagrangian:
            AdvectScalarsSemiLagrangian3D(*this, dt);
            break;
        case AdvectionScheme3D::MacCormackBFECC:
            AdvectScalarsMacCormackBFECC3D(*this, dt);
            break;
    }
    last_step_timings_.scalar_advection += SecondsSince(stage_start);

    stage_start = Clock::now();
    ApplyScalarAdvectionPostProcess_(dt);
    last_step_timings_.scalar_postprocess += SecondsSince(stage_start);
}

void EulerFluid3D::ApplyScalarAdvectionPostProcess_(Real dt) {
    const Real smoke_decay = std::exp(-std::max(0.0, config_.smoke_dissipation) * dt);
    const Real temp_decay = std::exp(-std::max(0.0, config_.temperature_dissipation) * dt);

    for (int x = 0; x < config_.resolution_x; ++x) {
        for (int y = 0; y < config_.resolution_y; ++y) {
            for (int z = 0; z < config_.resolution_z; ++z) {
                smoke_density_tmp_(x, y, z) = std::clamp(
                    smoke_density_tmp_(x, y, z) * smoke_decay,
                    0.0,
                    config_.max_smoke_density);
                temperature_tmp_(x, y, z) =
                    config_.ambient_temperature +
                    (temperature_tmp_(x, y, z) - config_.ambient_temperature) * temp_decay;
            }
        }
    }

    smoke_density_.Data() = smoke_density_tmp_.Data();
    temperature_.Data() = temperature_tmp_.Data();
}

void EulerFluid3D::ApplyConfiguredScalarSource_(Real dt) {
    if (!config_.source_enabled) {
        return;
    }

    AddSmokeSourceNormalized_(
        config_.source_min_x,
        config_.source_max_x,
        config_.source_min_y,
        config_.source_max_y,
        config_.source_min_z,
        config_.source_max_z,
        config_.source_smoke_rate,
        config_.source_temperature,
        config_.source_temperature_rate,
        dt);
}

void EulerFluid3D::ApplyConfiguredVelocitySource_(Real dt) {
    if (!config_.source_enabled) {
        return;
    }

    AddVelocityImpulseNormalized_(
        config_.source_min_x,
        config_.source_max_x,
        config_.source_min_y,
        config_.source_max_y,
        config_.source_min_z,
        config_.source_max_z,
        {config_.source_acceleration_x, config_.source_acceleration_y, config_.source_acceleration_z},
        dt);
}

void EulerFluid3D::ProjectVelocity_(Real dt) {
    if (dt <= 0.0) {
        return;
    }

    auto stage_start = Clock::now();
    ComputeDivergence_();
    last_step_timings_.projection_divergence += SecondsSince(stage_start);

    stage_start = Clock::now();
    const int width = config_.resolution_x;
    const int height = config_.resolution_y;
    const int depth = config_.resolution_z;
    const int cell_count = width * height * depth;
    const Real scale = dt / (config_.fluid_density * config_.cell_size * config_.cell_size);

    const Real rhs_sum = std::accumulate(
        divergence_.Data().begin(),
        divergence_.Data().end(),
        0.0,
        [](Real sum, Real div) {
            return sum - div;
        });
    const Real rhs_mean = rhs_sum / static_cast<Real>(cell_count);

    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(static_cast<std::size_t>(cell_count) * 7);
    Eigen::VectorXd rhs(cell_count);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                const int row = FlattenCellIndex3D(x, y, z, height, depth);
                if (row == 0) {
                    triplets.emplace_back(row, row, 1.0);
                    rhs(row) = 0.0;
                    continue;
                }

                rhs(row) = -divergence_(x, y, z) - rhs_mean;
                Real diagonal = 0.0;
                const int offsets[6][3] = {
                    {-1, 0, 0},
                    {1, 0, 0},
                    {0, -1, 0},
                    {0, 1, 0},
                    {0, 0, -1},
                    {0, 0, 1},
                };

                for (const auto& offset : offsets) {
                    const int nx = x + offset[0];
                    const int ny = y + offset[1];
                    const int nz = z + offset[2];
                    if (nx < 0 || nx >= width || ny < 0 || ny >= height || nz < 0 || nz >= depth) {
                        continue;
                    }

                    diagonal += scale;
                    const int neighbor = FlattenCellIndex3D(nx, ny, nz, height, depth);
                    if (neighbor != 0) {
                        triplets.emplace_back(row, neighbor, -scale);
                    }
                }

                triplets.emplace_back(row, row, diagonal > 0.0 ? diagonal : 1.0);
            }
        }
    }

    Eigen::SparseMatrix<Real> matrix(cell_count, cell_count);
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    last_step_timings_.pressure_system_build += SecondsSince(stage_start);

    stage_start = Clock::now();
    const PressureSolveResult3D pressure_result = SolvePressurePoisson3D(
        config_.pressure_solver,
        matrix,
        rhs,
        config_.pressure_iterations,
        config_.pressure_tolerance);
    last_step_timings_.pressure_solve += SecondsSince(stage_start);
    if (!pressure_result.success) {
        last_pressure_iterations_ = pressure_result.iterations;
        last_pressure_residual_ = pressure_result.residual;
        pressure_.Fill(0.0);
        return;
    }

    last_pressure_iterations_ = pressure_result.iterations;
    last_pressure_residual_ = pressure_result.residual;

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                pressure_(x, y, z) = pressure_result.solution(FlattenCellIndex3D(x, y, z, height, depth));
            }
        }
    }

    const std::vector<Real> u_before_projection = velocity_.UData();
    const std::vector<Real> v_before_projection = velocity_.VData();
    const std::vector<Real> w_before_projection = velocity_.WData();
    const Real grad_scale = dt / (config_.fluid_density * config_.cell_size);

    stage_start = Clock::now();
    for (int x = 1; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                velocity_.U(x, y, z) -= grad_scale * (pressure_(x, y, z) - pressure_(x - 1, y, z));
            }
        }
    }
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            for (int z = 0; z < depth; ++z) {
                velocity_.V(x, y, z) -= grad_scale * (pressure_(x, y, z) - pressure_(x, y - 1, z));
            }
        }
    }
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            for (int z = 1; z < depth; ++z) {
                velocity_.W(x, y, z) -= grad_scale * (pressure_(x, y, z) - pressure_(x, y, z - 1));
            }
        }
    }
    last_step_timings_.pressure_apply_gradient += SecondsSince(stage_start);

    stage_start = Clock::now();
    SetBoundaryConditions_();
    const Real projection_boundary_seconds = SecondsSince(stage_start);
    last_step_timings_.projection_boundary += projection_boundary_seconds;
    last_step_timings_.boundary_conditions += projection_boundary_seconds;

    stage_start = Clock::now();
    ComputeDivergence_();
    last_step_timings_.projection_divergence += SecondsSince(stage_start);

    stage_start = Clock::now();
    const Real max_velocity_after_projection = ComputeMaxSpeed_();
    if (!std::isfinite(max_velocity_after_projection) || max_velocity_after_projection > 1.0e6) {
        velocity_.UData() = u_before_projection;
        velocity_.VData() = v_before_projection;
        velocity_.WData() = w_before_projection;
        pressure_.Fill(0.0);
        last_pressure_residual_ = std::numeric_limits<Real>::infinity();
        SetBoundaryConditions_();
        ComputeDivergence_();
    }
    last_step_timings_.projection_validation += SecondsSince(stage_start);
}

void EulerFluid3D::SetBoundaryConditions_() {
    for (int y = 0; y < velocity_.UHeight(); ++y) {
        for (int z = 0; z < velocity_.UDepth(); ++z) {
            velocity_.U(0, y, z) = 0.0;
            velocity_.U(velocity_.UWidth() - 1, y, z) = 0.0;
        }
    }
    for (int x = 0; x < velocity_.VWidth(); ++x) {
        for (int z = 0; z < velocity_.VDepth(); ++z) {
            velocity_.V(x, 0, z) = 0.0;
            velocity_.V(x, velocity_.VHeight() - 1, z) = 0.0;
        }
    }
    for (int x = 0; x < velocity_.WWidth(); ++x) {
        for (int y = 0; y < velocity_.WHeight(); ++y) {
            velocity_.W(x, y, 0) = 0.0;
            velocity_.W(x, y, velocity_.WDepth() - 1) = 0.0;
        }
    }
}

void EulerFluid3D::ComputeDivergence_() {
    for (int x = 0; x < config_.resolution_x; ++x) {
        for (int y = 0; y < config_.resolution_y; ++y) {
            for (int z = 0; z < config_.resolution_z; ++z) {
                divergence_(x, y, z) =
                    (velocity_.U(x + 1, y, z) - velocity_.U(x, y, z) +
                     velocity_.V(x, y + 1, z) - velocity_.V(x, y, z) +
                     velocity_.W(x, y, z + 1) - velocity_.W(x, y, z)) *
                    inverse_cell_size_;
            }
        }
    }
}

Vector3R EulerFluid3D::TraceAdvectionPosition_(const Vector3R& position, Real dt) const {
    const auto velocity_at = [&](const Vector3R& p) {
        return SampleVelocityWorld_(velocity_, p.x(), p.y(), p.z());
    };
    return ClampWorldPosition_(BacktracePosition3D(config_.advection_integrator, position, dt, velocity_at));
}

Vector3R EulerFluid3D::ClampWorldPosition_(const Vector3R& position) const {
    return {
        std::clamp(position.x(), 0.0, static_cast<Real>(config_.resolution_x) * config_.cell_size),
        std::clamp(position.y(), 0.0, static_cast<Real>(config_.resolution_y) * config_.cell_size),
        std::clamp(position.z(), 0.0, static_cast<Real>(config_.resolution_z) * config_.cell_size),
    };
}

Real EulerFluid3D::SampleScalarWorld_(const CellCenteredScalarGrid3D& grid, Real x, Real y, Real z) const {
    return grid.Sample(
        CellCenterGridCoordinate(x, inverse_cell_size_),
        CellCenterGridCoordinate(y, inverse_cell_size_),
        CellCenterGridCoordinate(z, inverse_cell_size_),
        config_.advection_interpolation);
}

Vector3R EulerFluid3D::SampleVelocityWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const {
    return {SampleUWorld_(grid, x, y, z), SampleVWorld_(grid, x, y, z), SampleWWorld_(grid, x, y, z)};
}

Real EulerFluid3D::SampleUWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const {
    return grid.SampleU(
        FaceGridCoordinate(x, inverse_cell_size_),
        CellCenterGridCoordinate(y, inverse_cell_size_),
        CellCenterGridCoordinate(z, inverse_cell_size_),
        config_.advection_interpolation);
}

Real EulerFluid3D::SampleVWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const {
    return grid.SampleV(
        CellCenterGridCoordinate(x, inverse_cell_size_),
        FaceGridCoordinate(y, inverse_cell_size_),
        CellCenterGridCoordinate(z, inverse_cell_size_),
        config_.advection_interpolation);
}

Real EulerFluid3D::SampleWWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const {
    return grid.SampleW(
        CellCenterGridCoordinate(x, inverse_cell_size_),
        CellCenterGridCoordinate(y, inverse_cell_size_),
        FaceGridCoordinate(z, inverse_cell_size_),
        config_.advection_interpolation);
}

void EulerFluid3D::AddSmokeSourceNormalized_(
    Real min_x,
    Real max_x,
    Real min_y,
    Real max_y,
    Real min_z,
    Real max_z,
    Real smoke_rate,
    Real target_temperature,
    Real temperature_rate,
    Real dt) {
    const int x_begin = static_cast<int>(std::floor(min_x * config_.resolution_x));
    const int x_end = static_cast<int>(std::ceil(max_x * config_.resolution_x));
    const int y_begin = static_cast<int>(std::floor(min_y * config_.resolution_y));
    const int y_end = static_cast<int>(std::ceil(max_y * config_.resolution_y));
    const int z_begin = static_cast<int>(std::floor(min_z * config_.resolution_z));
    const int z_end = static_cast<int>(std::ceil(max_z * config_.resolution_z));
    const Real heat_alpha = 1.0 - std::exp(-std::max(0.0, temperature_rate) * dt);

    for (int x = x_begin; x < x_end; ++x) {
        for (int y = y_begin; y < y_end; ++y) {
            for (int z = z_begin; z < z_end; ++z) {
                smoke_density_(x, y, z) = std::clamp(
                    smoke_density_(x, y, z) + smoke_rate * dt,
                    0.0,
                    config_.max_smoke_density);
                temperature_(x, y, z) += heat_alpha * (target_temperature - temperature_(x, y, z));
            }
        }
    }
}

void EulerFluid3D::AddVelocityImpulseNormalized_(
    Real min_x,
    Real max_x,
    Real min_y,
    Real max_y,
    Real min_z,
    Real max_z,
    Vector3R acceleration,
    Real dt) {
    const Real domain_x = static_cast<Real>(config_.resolution_x) * config_.cell_size;
    const Real domain_y = static_cast<Real>(config_.resolution_y) * config_.cell_size;
    const Real domain_z = static_cast<Real>(config_.resolution_z) * config_.cell_size;
    const Real velocity_delta_x = acceleration.x() * dt * domain_x;
    const Real velocity_delta_y = acceleration.y() * dt * domain_y;
    const Real velocity_delta_z = acceleration.z() * dt * domain_z;

    for (int x = 0; x < velocity_.UWidth(); ++x) {
        for (int y = 0; y < velocity_.UHeight(); ++y) {
            for (int z = 0; z < velocity_.UDepth(); ++z) {
                const Real px = static_cast<Real>(x) / static_cast<Real>(config_.resolution_x);
                const Real py = (static_cast<Real>(y) + 0.5) / static_cast<Real>(config_.resolution_y);
                const Real pz = (static_cast<Real>(z) + 0.5) / static_cast<Real>(config_.resolution_z);
                if (px >= min_x && px <= max_x && py >= min_y && py <= max_y && pz >= min_z && pz <= max_z) {
                    velocity_.U(x, y, z) += velocity_delta_x;
                }
            }
        }
    }
    for (int x = 0; x < velocity_.VWidth(); ++x) {
        for (int y = 0; y < velocity_.VHeight(); ++y) {
            for (int z = 0; z < velocity_.VDepth(); ++z) {
                const Real px = (static_cast<Real>(x) + 0.5) / static_cast<Real>(config_.resolution_x);
                const Real py = static_cast<Real>(y) / static_cast<Real>(config_.resolution_y);
                const Real pz = (static_cast<Real>(z) + 0.5) / static_cast<Real>(config_.resolution_z);
                if (px >= min_x && px <= max_x && py >= min_y && py <= max_y && pz >= min_z && pz <= max_z) {
                    velocity_.V(x, y, z) += velocity_delta_y;
                }
            }
        }
    }
    for (int x = 0; x < velocity_.WWidth(); ++x) {
        for (int y = 0; y < velocity_.WHeight(); ++y) {
            for (int z = 0; z < velocity_.WDepth(); ++z) {
                const Real px = (static_cast<Real>(x) + 0.5) / static_cast<Real>(config_.resolution_x);
                const Real py = (static_cast<Real>(y) + 0.5) / static_cast<Real>(config_.resolution_y);
                const Real pz = static_cast<Real>(z) / static_cast<Real>(config_.resolution_z);
                if (px >= min_x && px <= max_x && py >= min_y && py <= max_y && pz >= min_z && pz <= max_z) {
                    velocity_.W(x, y, z) += velocity_delta_z;
                }
            }
        }
    }
}

}  // namespace skyspaces::reference
