#include "AdvectionIntegrator3D.h"

namespace skyspaces::reference {

Vector3R BacktracePosition3D(
    AdvectionIntegrator3D integrator,
    const Vector3R& position,
    Real dt,
    const std::function<Vector3R(const Vector3R&)>& velocity_at) {
    switch (integrator) {
        case AdvectionIntegrator3D::ExplicitEuler: {
            const Vector3R k1 = velocity_at(position);
            return position - dt * k1;
        }
        case AdvectionIntegrator3D::Midpoint: {
            const Vector3R k1 = velocity_at(position);
            const Vector3R midpoint = position - 0.5 * dt * k1;
            const Vector3R k2 = velocity_at(midpoint);
            return position - dt * k2;
        }
        case AdvectionIntegrator3D::RungeKutta3: {
            const Vector3R k1 = -velocity_at(position);
            const Vector3R k2 = -velocity_at(position + 0.5 * dt * k1);
            const Vector3R k3 = -velocity_at(position + 0.75 * dt * k2);
            return position + dt * ((2.0 / 9.0) * k1 + (3.0 / 9.0) * k2 + (4.0 / 9.0) * k3);
        }
        case AdvectionIntegrator3D::RungeKutta4: {
            const Vector3R k1 = -velocity_at(position);
            const Vector3R k2 = -velocity_at(position + 0.5 * dt * k1);
            const Vector3R k3 = -velocity_at(position + 0.5 * dt * k2);
            const Vector3R k4 = -velocity_at(position + dt * k3);
            return position + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
    }

    const Vector3R k1 = velocity_at(position);
    return position - dt * k1;
}

}  // namespace skyspaces::reference
