#include "PressureSolver2D.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <algorithm>
#include <limits>

namespace skyspaces {

namespace {

Real RelativeResidual(
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& rhs) {
    const Real rhs_norm = rhs.norm();
    const Real denominator = std::max(rhs_norm, std::numeric_limits<Real>::epsilon());
    return (matrix * x - rhs).norm() / denominator;
}

template <typename IterativeSolver>
PressureSolveResult2D SolveIterative(
    IterativeSolver& solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs,
    int max_iterations,
    Real tolerance) {
    PressureSolveResult2D result;
    result.solution = Eigen::VectorXd::Zero(rhs.size());

    solver.setMaxIterations(std::max(1, max_iterations));
    solver.setTolerance(tolerance);
    solver.compute(matrix);
    if (solver.info() != Eigen::Success) {
        result.residual = std::numeric_limits<Real>::infinity();
        return result;
    }

    result.solution = solver.solve(rhs);
    result.iterations = static_cast<int>(solver.iterations());
    result.residual = RelativeResidual(matrix, result.solution, rhs);
    result.success = result.solution.allFinite() && result.residual <= tolerance;
    return result;
}

template <typename DirectSolver>
PressureSolveResult2D SolveDirect(
    DirectSolver& solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs) {
    PressureSolveResult2D result;
    result.solution = Eigen::VectorXd::Zero(rhs.size());

    solver.compute(matrix);
    if (solver.info() != Eigen::Success) {
        result.residual = std::numeric_limits<Real>::infinity();
        return result;
    }

    result.solution = solver.solve(rhs);
    result.iterations = 1;
    result.residual = RelativeResidual(matrix, result.solution, rhs);
    result.success = solver.info() == Eigen::Success;
    return result;
}

}  // namespace

PressureSolveResult2D SolvePressurePoisson2D(
    PressureSolver2D solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs,
    int max_iterations,
    Real tolerance) {
    switch (solver) {
        case PressureSolver2D::ConjugateGradient: {
            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<Real>,
                Eigen::Lower | Eigen::Upper> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::BiCGSTAB: {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::SimplicialLDLT: {
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveDirect(eigen_solver, matrix, rhs);
        }
        case PressureSolver2D::SparseLU: {
            Eigen::SparseLU<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveDirect(eigen_solver, matrix, rhs);
        }
    }

    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<Real>,
        Eigen::Lower | Eigen::Upper> eigen_solver;
    return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
}

}  // namespace skyspaces
