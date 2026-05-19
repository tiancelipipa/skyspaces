#pragma once

#include <vector>
#include "Vector2D.h"
#include "MACGrid2D.h"

namespace skyspaces {

class Fluid2D {
public:
    struct Config {
        int resolution_x = 128;
        int resolution_y = 128;
        
        double cell_size = 1.0;
        double time_step = 1e-3;
        double cfl_number = 0.75;   // time_step will be automatically computed based on this and max speed if time_step <= 0
        int max_substeps = 8;
        int semi_lagrangian_substeps = 1;
        
        
        // double pressure_tolerance = 1e-5;
        // int pressure_iterations = 2000;
        // int diffusion_iterations = 40;
        // double density_dissipation = 1.0;
        // double velocity_dissipation = 1.0;
        // double viscosity = 5e-4;
        // double buoyancy = 0.1;
        // double vorticity_confinement = 0.12;
        // bool wind_enabled = false;
        // double wind_strength = 2e-2;
        

        
        // bool source_enabled = true;
        // double source_min_x = 0.45;
        // double source_max_x = 0.55;
        // double source_min_y = 0.10;
        // double source_max_y = 0.15;
        // double source_density = 1.0;
    };

    explicit Fluid2D(Config config = {});

    void Reset();
    void Step();
    void Step(double dt);

    // void AddDensitySource(
    //     double min_x,
    //     double max_x,
    //     double min_y,
    //     double max_y,
    //     double density = 1.0);

    // void AddForce(
    //     double min_x,
    //     double max_x,
    //     double min_y,
    //     double max_y,
    //     double force_x,
    //     double force_y);

    // void ClearDensity();
    // void ClearForces();

    // int ResolutionX() const noexcept;
    // int ResolutionY() const noexcept;
    // double CellSize() const noexcept;
    // double InverseCellSize() const noexcept;
    // double DomainSizeX() const noexcept;
    // double DomainSizeY() const noexcept;
    // double TimeStep() const noexcept;
    // double Time() const noexcept;
    // int StepCount() const noexcept;
    // int LastPressureIterations() const noexcept;
    // double LastPressureResidual() const noexcept;

    // Vector2D SampleVelocityNormalized(double x, double y) const;
    // double SampleDensityNormalized(double x, double y) const;

    // const Config& Configs() const noexcept;
    // Config& Configs() noexcept;

    // const MACScalarGrid2D& Density() const noexcept;
    // const MACScalarGrid2D& Pressure() const noexcept;
    // const MACScalarGrid2D& Divergence() const noexcept;
    // const MACScalarGrid2D& Vorticity() const noexcept;
    // const MACScalarGrid2D& VelocityX() const noexcept;
    // const MACScalarGrid2D& VelocityY() const noexcept;

private:
    void InitializeFromConfig_();
    int ComputeSubstepCount_(double dt) const;
    double ComputeMaxSpeed_() const;

    void AdvectFields_(double dt);
    void SemiLagrangianBackwards_(double dt, int substeps = 1);
    // void ApplyConfiguredSource();
    // void ApplyBuoyancy();
    // void ApplyWind();
    // void ApplyVorticityConfinement();
    // void ApplyAccumulatedForces(double dt);
    // void DiffuseVelocity(double dt);

    // void SolvePressureProjection(double dt);
    // void SetBoundaryConditions();
    // void ComputeDivergenceField();
    // void SolvePoisson(double dt);
    // void CorrectVelocity(double dt);
    // void ComputeVorticityField();
    
    

    // double SampleDensityAtCenterCoords(double x, double y) const;
    // double SampleUAtCenterCoords(double x, double y) const;
    // double SampleVAtCenterCoords(double x, double y) const;
    // Vector2D SampleVelocityAtCenterCoords(double x, double y) const;
    // Vector2D SampleVelocityAtUFaceCoords(double x, double y) const;
    // Vector2D SampleVelocityAtVFaceCoords(double x, double y) const;
    // double PressureAtWithNeumann(int x, int y) const;

    Config config_;

    double inverse_cell_size_ = 1.0;
    double time_ = 0.0;
    int step_count_ = 0;
    // double domain_size_y_ = 1.0;
    // double wind_phase_ = 0.0;
    
    // int last_pressure_iterations_ = 0;
    // double last_pressure_residual_ = 0.0;

    // MACScalarGrid2D density_;
    // MACScalarGrid2D density_tmp_;
    // MACScalarGrid2D pressure_;
    // MACScalarGrid2D divergence_;
    // MACScalarGrid2D vorticity_;
    MACVectorGrid2D velocity_;
    // MACScalarGrid2D velocity_x_tmp_;
    // MACScalarGrid2D velocity_y_tmp_;
    MACVectorGrid2D semi_lagrangian_position_;
    MACVectorGrid2D force_;
    // MACScalarGrid2D external_force_x_;
    // MACScalarGrid2D external_force_y_;
};

}  // namespace skyspaces
