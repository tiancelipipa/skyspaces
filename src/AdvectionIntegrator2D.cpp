#include "AdvectionIntegrator2D.h"

namespace skyspaces {

Vector2D BacktracePosition2D(
    AdvectionIntegrator2D integrator,
    const Vector2D& position,
    Real dt,
    const std::function<Vector2D(const Vector2D&)>& velocity_at) {
    switch (integrator) {
        case AdvectionIntegrator2D::ExplicitEuler: {
            const Vector2D k1 = velocity_at(position);
            return position - dt * k1;
        }
        case AdvectionIntegrator2D::Midpoint: {
            const Vector2D k1 = velocity_at(position);
            const Vector2D midpoint = position - 0.5 * dt * k1;
            const Vector2D k2 = velocity_at(midpoint);
            return position - dt * k2;
        }
        case AdvectionIntegrator2D::RungeKutta3: {
            const Vector2D k1 = -velocity_at(position);
            const Vector2D k2 = -velocity_at(position + 0.5 * dt * k1);
            const Vector2D k3 = -velocity_at(position + 0.75 * dt * k2);
            return position + dt * ((2.0 / 9.0) * k1 + (3.0 / 9.0) * k2 + (4.0 / 9.0) * k3);
        }
        case AdvectionIntegrator2D::RungeKutta4: {
            const Vector2D k1 = -velocity_at(position);
            const Vector2D k2 = -velocity_at(position + 0.5 * dt * k1);
            const Vector2D k3 = -velocity_at(position + 0.5 * dt * k2);
            const Vector2D k4 = -velocity_at(position + dt * k3);
            return position + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
    }

    const Vector2D k1 = velocity_at(position);
    return position - dt * k1;
}

}  // namespace skyspaces
