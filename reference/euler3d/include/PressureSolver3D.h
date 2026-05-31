#pragma once

#include "MathTypes.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace skyspaces::reference {

enum class PressureSolver3D {
    Jacobi,
    GaussSeidel,
    SuccessiveOverRelaxation,
    ConjugateGradient,
    ConjugateGradientNoPreconditioner,
    BiCGSTAB,
    LeastSquaresConjugateGradient,
    SimplicialLLT,
    SimplicialLDLT,
    SparseLU,
};

struct PressureSolveResult3D {
    Eigen::VectorXd solution;
    int iterations = 0;
    Real residual = 0.0;
    bool success = false;
};

PressureSolveResult3D SolvePressurePoisson3D(
    PressureSolver3D solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs,
    int max_iterations,
    Real tolerance);

}  // namespace skyspaces::reference
