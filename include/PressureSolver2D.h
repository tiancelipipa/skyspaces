#pragma once

#include "MathTypes.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace skyspaces {

enum class PressureSolver2D {
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

struct PressureSolveResult2D {
    Eigen::VectorXd solution;
    int iterations = 0;
    Real residual = 0.0;
    bool success = false;
};

PressureSolveResult2D SolvePressurePoisson2D(
    PressureSolver2D solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs,
    int max_iterations,
    Real tolerance);

}  // namespace skyspaces
