#pragma once

#include "MathTypes.h"

#include <functional>

namespace skyspaces {

enum class AdvectionIntegrator2D {
    ExplicitEuler,
    Midpoint,
    RungeKutta3,
    RungeKutta4,
};

Vector2R BacktracePosition2D(
    AdvectionIntegrator2D integrator,
    const Vector2R& position,
    Real dt,
    const std::function<Vector2R(const Vector2R&)>& velocity_at);

}  // namespace skyspaces
