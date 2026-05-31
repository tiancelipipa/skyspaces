#pragma once

#include "Types3D.h"

#include <functional>

namespace skyspaces::reference {

enum class AdvectionIntegrator3D {
    ExplicitEuler,
    Midpoint,
    RungeKutta3,
    RungeKutta4,
};

Vector3R BacktracePosition3D(
    AdvectionIntegrator3D integrator,
    const Vector3R& position,
    Real dt,
    const std::function<Vector3R(const Vector3R&)>& velocity_at);

}  // namespace skyspaces::reference
