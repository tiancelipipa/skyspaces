#pragma once

#include "AdvectionIntegrator2D.h"
#include "AdvectionScheme2D.h"
#include "Grid2D.h"
#include "PressureSolver2D.h"
#include "Solid2D.h"

#include <optional>

namespace skyspaces {

struct Fluid2DAdvection2D;

// Small 2D smoke simulation using a staggered velocity grid and scalar
// cell-centered fields. The solver favors clarity and stable building blocks
// over high-order advection or production renderer features.
class Fluid2D {
public:
    using SolidBoundaryMode = SolidBoundaryMode2D;

    friend struct Fluid2DAdvection2D;

    struct Config {
        // Small positive threshold used to avoid division by near-zero values.
        Real numeric_epsilon = 1e-12;
        
        int resolution_x = 128;
        int resolution_y = 128;
        Real cell_size = 1.0;

        Real time_step = 1.0 / 60.0;
        int max_substeps = 8;
        // CFL([1.0,10.0]) substepping keeps semi-Lagrangian backtracing
        // from taking very large jumps when velocities spike.
        Real cfl_number = 1.0;

        AdvectionIntegrator2D advection_integrator = AdvectionIntegrator2D::ExplicitEuler;
        InterpolationMethod2D advection_interpolation = InterpolationMethod2D::Bilinear;
        AdvectionScheme2D advection_scheme = AdvectionScheme2D::SemiLagrangian;
        PressureSolver2D pressure_solver = PressureSolver2D::ConjugateGradient;

        Real fluid_density = 1.0;           // 流体介质的物理密度，影响压力/速度求解
        int pressure_iterations = 400;
        Real pressure_tolerance = 1e-6;

        Real ambient_temperature = 273.15;
        // Vorticity confinement is optional; it visually restores small-scale
        // swirl lost to advection.
        Real vorticity_confinement = 0.0;

        Real smoke_dissipation = 0.05;
        Real temperature_dissipation = 0.10;
        Real max_smoke_density = 1.0;       // 烟雾浓度，影响画面里烟有多浓，与 fluid_density 不同

        // Source bounds are normalized to the simulation domain [0, 1]^2.
        bool source_enabled = true;
        Real source_min_x = 0.45;
        Real source_max_x = 0.55;
        Real source_min_y = 0.05;
        Real source_max_y = 0.12;
        Real source_smoke_rate = 2.0;
        Real source_temperature = 650.0;
        Real source_temperature_rate = 6.0;
        Real source_acceleration_x = 1.5;   // domain widths per second^2
        Real source_acceleration_y = 0.0;   // domain heights per second^2
    };

    explicit Fluid2D(Config config = {});

    void Reset();
    void Step();
    void Step(Real dt);

    void ClearSolidBoundary();
    void SetSolidCellMarkers(const CellCenteredScalarGrid2D& solid_cells);
    void SetSolidLevelSet(const CellCenteredScalarGrid2D& solid_phi);
    void SetSolidVelocity(const FaceCenteredVectorGrid2D& solid_velocity);

    int ResolutionX() const noexcept;
    int ResolutionY() const noexcept;
    Real CellSize() const noexcept;
    Real TimeStep() const noexcept;
    Real Time() const noexcept;
    int StepCount() const noexcept;
    int LastPressureIterations() const noexcept;
    Real LastPressureResidual() const noexcept;

    // Sampling helpers accept normalized coordinates in [0, 1]^2.
    Real SampleSmokeNormalized(Real x, Real y) const;
    Real SampleTemperatureNormalized(Real x, Real y) const;
    Vector2R SampleVelocityNormalized(Real x, Real y) const;

    const Config& Configs() const noexcept;
    Config& Configs() noexcept;

    const CellCenteredScalarGrid2D& SmokeDensity() const noexcept;
    CellCenteredScalarGrid2D& SmokeDensity() noexcept;
    const CellCenteredScalarGrid2D& Temperature() const noexcept;
    CellCenteredScalarGrid2D& Temperature() noexcept;
    const CellCenteredScalarGrid2D& Pressure() const noexcept;
    CellCenteredScalarGrid2D& Pressure() noexcept;
    const CellCenteredScalarGrid2D& Divergence() const noexcept;
    const CellCenteredScalarGrid2D& Vorticity() const noexcept;
    const CellCenteredScalarGrid2D& SolidCellMarkers() const noexcept;
    const CellCenteredScalarGrid2D& SolidLevelSet() const noexcept;
    const FaceCenteredVectorGrid2D& Velocity() const noexcept;
    const FaceCenteredVectorGrid2D& SolidVelocity() const noexcept;
    const Solid2D& Solid() const noexcept;

private:
    void InitializeFromConfig_();
    void ConfigureAdvectionScratch_();
    void ClearAdvectionScratch_();
    int ComputeSubstepCount_(Real dt) const;
    Real ComputeMaxSpeed_() const;

    void ApplyConfiguredScalarSource_(Real dt);
    void ApplyConfiguredVelocitySource_(Real dt);
    void ComputeVorticity_();
    void ApplyVorticityConfinement_(Real dt);

    void AdvectVelocity_(Real dt);
    void AdvectScalars_(Real dt);
    void ApplyScalarAdvectionPostProcess_(Real dt);
    void ProjectVelocity_(Real dt);
    void SetBoundaryConditions_();
    void ComputeDivergence_();

    bool HasSolidBoundary_() const noexcept;
    bool IsSolidCell_(int x, int y) const;
    bool IsSolidWorld_(Real x, Real y) const;
    bool IsUFaceOpen_(int x, int y) const;
    bool IsVFaceOpen_(int x, int y) const;
    Real SolidU_(int x, int y) const;
    Real SolidV_(int x, int y) const;
    Vector2R ProjectOutOfSolid_(const Vector2R& position, const Vector2R& fallback) const;

    Real SampleScalarWorld_(const CellCenteredScalarGrid2D& grid, Real x, Real y) const;
    Vector2R SampleVelocityWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const;
    Real SampleUWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const;
    Real SampleVWorld_(const FaceCenteredVectorGrid2D& grid, Real x, Real y) const;

    // Adds smoke and heat inside a normalized axis-aligned rectangle.
    void AddSmokeSourceNormalized_(
        Real min_x,
        Real max_x,
        Real min_y,
        Real max_y,
        Real smoke_rate,
        Real target_temperature,
        Real temperature_rate,
        Real dt);

    // Adds a velocity delta on staggered faces overlapped by the normalized box.
    void AddVelocityImpulseNormalized_(
        Real min_x,
        Real max_x,
        Real min_y,
        Real max_y,
        Vector2R acceleration,
        Real dt);

    Config config_;

    Real inverse_cell_size_ = 1.0;
    Real time_ = 0.0;
    int step_count_ = 0;
    int last_pressure_iterations_ = 0;
    Real last_pressure_residual_ = 0.0;

    CellCenteredScalarGrid2D smoke_density_;
    CellCenteredScalarGrid2D smoke_density_tmp_;
    CellCenteredScalarGrid2D temperature_;
    CellCenteredScalarGrid2D temperature_tmp_;
    CellCenteredScalarGrid2D pressure_;
    CellCenteredScalarGrid2D divergence_;
    CellCenteredScalarGrid2D vorticity_;
    FaceCenteredVectorGrid2D velocity_;
    FaceCenteredVectorGrid2D velocity_tmp_;

    std::optional<FaceCenteredVectorGrid2D> velocity_first_pass_;
    std::optional<FaceCenteredVectorGrid2D> velocity_back_pass_;
    std::optional<FaceCenteredVectorGrid2D> velocity_corrected_source_;
    std::optional<CellCenteredScalarGrid2D> smoke_first_pass_;
    std::optional<CellCenteredScalarGrid2D> smoke_back_pass_;
    std::optional<CellCenteredScalarGrid2D> smoke_corrected_source_;
    std::optional<CellCenteredScalarGrid2D> temperature_first_pass_;
    std::optional<CellCenteredScalarGrid2D> temperature_back_pass_;
    std::optional<CellCenteredScalarGrid2D> temperature_corrected_source_;

    Solid2D solid_;
};

}  // namespace skyspaces
