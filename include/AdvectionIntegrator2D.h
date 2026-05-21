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

Vector2D BacktracePosition2D(
    AdvectionIntegrator2D integrator,
    const Vector2D& position,
    Real dt,
    const std::function<Vector2D(const Vector2D&)>& velocity_at);

}  // namespace skyspaces
