#pragma once

#include "MathTypes.h"

namespace skyspaces {

class Fluid2D;

enum class AdvectionScheme2D {
    SemiLagrangian,
    MacCormackBFECC,
};

void AdvectVelocitySemiLagrangian2D(Fluid2D& fluid, Real dt);
void AdvectScalarsSemiLagrangian2D(Fluid2D& fluid, Real dt);
void AdvectVelocityMacCormackBFECC2D(Fluid2D& fluid, Real dt);
void AdvectScalarsMacCormackBFECC2D(Fluid2D& fluid, Real dt);

}  // namespace skyspaces
