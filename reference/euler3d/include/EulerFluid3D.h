#pragma once

#include "AdvectionIntegrator3D.h"
#include "AdvectionScheme3D.h"
#include "Interpolation3D.h"
#include "PressureSolver3D.h"
#include "Types3D.h"

#include <optional>
#include <vector>

namespace skyspaces::reference {

struct EulerFluid3DAdvection3D;

class CellCenteredScalarGrid3D {
public:
    CellCenteredScalarGrid3D() = default;
    CellCenteredScalarGrid3D(int width, int height, int depth, Real value = 0.0);

    void Resize(int width, int height, int depth, Real value = 0.0);
    void Fill(Real value);

    int Width() const noexcept;
    int Height() const noexcept;
    int Depth() const noexcept;
    bool Empty() const noexcept;

    Real& operator()(int x, int y, int z);
    Real operator()(int x, int y, int z) const;

    Real Sample(Real x, Real y, Real z) const;
    Real Sample(Real x, Real y, Real z, InterpolationMethod3D method) const;

    const std::vector<Real>& Data() const noexcept;
    std::vector<Real>& Data() noexcept;

private:
    int Index_(int x, int y, int z) const noexcept;

    int width_ = 0;
    int height_ = 0;
    int depth_ = 0;
    std::vector<Real> data_;
};

class FaceCenteredVectorGrid3D {
public:
    FaceCenteredVectorGrid3D() = default;
    FaceCenteredVectorGrid3D(int width, int height, int depth, Real value = 0.0);

    void Resize(int width, int height, int depth, Real value = 0.0);
    void Fill(Real value);

    int Width() const noexcept;
    int Height() const noexcept;
    int Depth() const noexcept;

    int UWidth() const noexcept;
    int UHeight() const noexcept;
    int UDepth() const noexcept;
    int VWidth() const noexcept;
    int VHeight() const noexcept;
    int VDepth() const noexcept;
    int WWidth() const noexcept;
    int WHeight() const noexcept;
    int WDepth() const noexcept;

    Real& U(int x, int y, int z);
    Real U(int x, int y, int z) const;
    Real& V(int x, int y, int z);
    Real V(int x, int y, int z) const;
    Real& W(int x, int y, int z);
    Real W(int x, int y, int z) const;

    Real SampleU(Real x, Real y, Real z) const;
    Real SampleU(Real x, Real y, Real z, InterpolationMethod3D method) const;
    Real SampleV(Real x, Real y, Real z) const;
    Real SampleV(Real x, Real y, Real z, InterpolationMethod3D method) const;
    Real SampleW(Real x, Real y, Real z) const;
    Real SampleW(Real x, Real y, Real z, InterpolationMethod3D method) const;
    Vector3R Sample(Real x, Real y, Real z) const;
    Vector3R Sample(Real x, Real y, Real z, InterpolationMethod3D method) const;

    const std::vector<Real>& UData() const noexcept;
    std::vector<Real>& UData() noexcept;
    const std::vector<Real>& VData() const noexcept;
    std::vector<Real>& VData() noexcept;
    const std::vector<Real>& WData() const noexcept;
    std::vector<Real>& WData() noexcept;

private:
    int UIndex_(int x, int y, int z) const noexcept;
    int VIndex_(int x, int y, int z) const noexcept;
    int WIndex_(int x, int y, int z) const noexcept;

    int width_ = 0;
    int height_ = 0;
    int depth_ = 0;
    std::vector<Real> u_data_;
    std::vector<Real> v_data_;
    std::vector<Real> w_data_;
};

class EulerFluid3D {
public:
    friend struct EulerFluid3DAdvection3D;

    struct StepTimings {
        int substep_count = 0;
        Real total = 0.0;
        Real advect_velocity = 0.0;
        Real velocity_source = 0.0;
        Real boundary_conditions = 0.0;
        Real pre_projection_boundary = 0.0;
        Real projection_boundary = 0.0;
        Real projection = 0.0;
        Real projection_divergence = 0.0;
        Real pressure_system_build = 0.0;
        Real pressure_solve = 0.0;
        Real pressure_apply_gradient = 0.0;
        Real projection_validation = 0.0;
        Real advect_scalars = 0.0;
        Real scalar_advection = 0.0;
        Real scalar_postprocess = 0.0;
        Real scalar_source = 0.0;
    };

    struct Config {
        Real numeric_epsilon = 1e-12;

        int resolution_x = 32;
        int resolution_y = 32;
        int resolution_z = 32;
        Real cell_size = 1.0 / 32.0;

        Real time_step = 1.0 / 60.0;
        int max_substeps = 8;
        Real cfl_number = 1.0;
        AdvectionIntegrator3D advection_integrator = AdvectionIntegrator3D::ExplicitEuler;
        InterpolationMethod3D advection_interpolation = InterpolationMethod3D::Trilinear;
        AdvectionScheme3D advection_scheme = AdvectionScheme3D::SemiLagrangian;
        PressureSolver3D pressure_solver = PressureSolver3D::ConjugateGradient;

        Real fluid_density = 1.0;
        int pressure_iterations = 400;
        Real pressure_tolerance = 1e-6;

        Real ambient_temperature = 273.15;
        Real smoke_dissipation = 0.05;
        Real temperature_dissipation = 0.10;
        Real max_smoke_density = 1.0;

        bool source_enabled = true;
        Real source_min_x = 0.45;
        Real source_max_x = 0.55;
        Real source_min_y = 0.05;
        Real source_max_y = 0.12;
        Real source_min_z = 0.45;
        Real source_max_z = 0.55;
        Real source_smoke_rate = 2.0;
        Real source_temperature = 650.0;
        Real source_temperature_rate = 6.0;
        Real source_acceleration_x = 0.0;
        Real source_acceleration_y = 2.0;
        Real source_acceleration_z = 0.0;
    };

    explicit EulerFluid3D(Config config = {});

    void Reset();
    void Step();
    void Step(Real dt);

    int ResolutionX() const noexcept;
    int ResolutionY() const noexcept;
    int ResolutionZ() const noexcept;
    Real CellSize() const noexcept;
    Real TimeStep() const noexcept;
    Real Time() const noexcept;
    int StepCount() const noexcept;
    int LastPressureIterations() const noexcept;
    Real LastPressureResidual() const noexcept;
    const StepTimings& LastStepTimings() const noexcept;

    Real SampleSmokeNormalized(Real x, Real y, Real z) const;
    Real SampleTemperatureNormalized(Real x, Real y, Real z) const;
    Vector3R SampleVelocityNormalized(Real x, Real y, Real z) const;

    const Config& Configs() const noexcept;
    Config& Configs() noexcept;

    const CellCenteredScalarGrid3D& SmokeDensity() const noexcept;
    CellCenteredScalarGrid3D& SmokeDensity() noexcept;
    const CellCenteredScalarGrid3D& Temperature() const noexcept;
    CellCenteredScalarGrid3D& Temperature() noexcept;
    const CellCenteredScalarGrid3D& Pressure() const noexcept;
    const CellCenteredScalarGrid3D& Divergence() const noexcept;
    const FaceCenteredVectorGrid3D& Velocity() const noexcept;

private:
    void InitializeFromConfig_();
    void ConfigureAdvectionScratch_();
    void ClearAdvectionScratch_();
    int ComputeSubstepCount_(Real dt) const;
    Real ComputeMaxSpeed_() const;

    void AdvectVelocity_(Real dt);
    void AdvectScalars_(Real dt);
    void ApplyScalarAdvectionPostProcess_(Real dt);
    void ApplyConfiguredScalarSource_(Real dt);
    void ApplyConfiguredVelocitySource_(Real dt);
    void ProjectVelocity_(Real dt);
    void SetBoundaryConditions_();
    void ComputeDivergence_();

    Vector3R TraceAdvectionPosition_(const Vector3R& position, Real dt) const;
    Vector3R ClampWorldPosition_(const Vector3R& position) const;
    Real SampleScalarWorld_(const CellCenteredScalarGrid3D& grid, Real x, Real y, Real z) const;
    Vector3R SampleVelocityWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const;
    Real SampleUWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const;
    Real SampleVWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const;
    Real SampleWWorld_(const FaceCenteredVectorGrid3D& grid, Real x, Real y, Real z) const;

    void AddSmokeSourceNormalized_(
        Real min_x,
        Real max_x,
        Real min_y,
        Real max_y,
        Real min_z,
        Real max_z,
        Real smoke_rate,
        Real target_temperature,
        Real temperature_rate,
        Real dt);

    void AddVelocityImpulseNormalized_(
        Real min_x,
        Real max_x,
        Real min_y,
        Real max_y,
        Real min_z,
        Real max_z,
        Vector3R acceleration,
        Real dt);

    Config config_;
    Real inverse_cell_size_ = 1.0;
    Real time_ = 0.0;
    int step_count_ = 0;
    int last_pressure_iterations_ = 0;
    Real last_pressure_residual_ = 0.0;
    StepTimings last_step_timings_;

    CellCenteredScalarGrid3D smoke_density_;
    CellCenteredScalarGrid3D smoke_density_tmp_;
    CellCenteredScalarGrid3D temperature_;
    CellCenteredScalarGrid3D temperature_tmp_;
    CellCenteredScalarGrid3D pressure_;
    CellCenteredScalarGrid3D divergence_;
    FaceCenteredVectorGrid3D velocity_;
    FaceCenteredVectorGrid3D velocity_tmp_;

    std::optional<FaceCenteredVectorGrid3D> velocity_first_pass_;
    std::optional<FaceCenteredVectorGrid3D> velocity_back_pass_;
    std::optional<FaceCenteredVectorGrid3D> velocity_corrected_source_;
    std::optional<CellCenteredScalarGrid3D> smoke_first_pass_;
    std::optional<CellCenteredScalarGrid3D> smoke_back_pass_;
    std::optional<CellCenteredScalarGrid3D> smoke_corrected_source_;
    std::optional<CellCenteredScalarGrid3D> temperature_first_pass_;
    std::optional<CellCenteredScalarGrid3D> temperature_back_pass_;
    std::optional<CellCenteredScalarGrid3D> temperature_corrected_source_;
};

}  // namespace skyspaces::reference
