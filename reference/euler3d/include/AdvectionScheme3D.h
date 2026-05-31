#pragma once

#include "MathTypes.h"

namespace skyspaces::reference {

class EulerFluid3D;

enum class AdvectionScheme3D {
    SemiLagrangian,
    MacCormackBFECC,
};

void AdvectVelocitySemiLagrangian3D(EulerFluid3D& fluid, Real dt);
void AdvectScalarsSemiLagrangian3D(EulerFluid3D& fluid, Real dt);
void AdvectVelocityMacCormackBFECC3D(EulerFluid3D& fluid, Real dt);
void AdvectScalarsMacCormackBFECC3D(EulerFluid3D& fluid, Real dt);

}  // namespace skyspaces::reference
