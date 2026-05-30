#include "PressureSolver2D.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>

namespace skyspaces {

namespace {

constexpr Real kSorOmega = 1.5;

Real RelativeResidual(
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& rhs) {
    const Real rhs_norm = rhs.norm();
    const Real denominator = std::max(rhs_norm, std::numeric_limits<Real>::epsilon());
    return (matrix * x - rhs).norm() / denominator;
}

enum class StationaryPressureSolver2D {
    Jacobi,
    GaussSeidel,
    SuccessiveOverRelaxation,
};

PressureSolveResult2D SolveStationary(
    StationaryPressureSolver2D solver,
    const Eigen::SparseMatrix<Real>& matrix,
    const Eigen::VectorXd& rhs,
    int max_iterations,
    Real tolerance) {
    using RowMajorSparseMatrix = Eigen::SparseMatrix<Real, Eigen::RowMajor>;

    const RowMajorSparseMatrix row_matrix = matrix;
    const int size = static_cast<int>(rhs.size());
    const int iteration_count = std::max(1, max_iterations);
    const Real diagonal_epsilon = std::numeric_limits<Real>::epsilon();

    PressureSolveResult2D result;
    result.solution = Eigen::VectorXd::Zero(size);
    Eigen::VectorXd next_solution = Eigen::VectorXd::Zero(size);

    for (int iteration = 0; iteration < iteration_count; ++iteration) {
        if (solver == StationaryPressureSolver2D::Jacobi) {
            for (int row = 0; row < size; ++row) {
                Real diagonal = 0.0;
                Real sum = rhs(row);
                for (RowMajorSparseMatrix::InnerIterator entry(row_matrix, row); entry; ++entry) {
                    if (entry.col() == row) {
                        diagonal = entry.value();
                    } else {
                        sum -= entry.value() * result.solution(entry.col());
                    }
                }

                next_solution(row) = std::abs(diagonal) > diagonal_epsilon
                    ? sum / diagonal
                    : result.solution(row);
            }
            result.solution.swap(next_solution);
        } else {
            const Real omega = solver == StationaryPressureSolver2D::SuccessiveOverRelaxation
                ? kSorOmega
                : 1.0;
            for (int row = 0; row < size; ++row) {
                Real diagonal = 0.0;
                Real sum = rhs(row);
                for (RowMajorSparseMatrix::InnerIterator entry(row_matrix, row); entry; ++entry) {
                    if (entry.col() == row) {
                        diagonal = entry.value();
                    } else {
                        sum -= entry.value() * result.solution(entry.col());
                    }
                }

                if (std::abs(diagonal) > diagonal_epsilon) {
                    const Real updated = sum / diagonal;
                    result.solution(row) += omega * (updated - result.solution(row));
                }
            }
        }

        result.iterations = iteration + 1;
        result.residual = RelativeResidual(matrix, result.solution, rhs);
        result.success =
            result.solution.allFinite() &&
            result.residual <= tolerance;
        if (result.success) {
            return result;
        }
    }

    result.residual = RelativeResidual(matrix, result.solution, rhs);
    result.success =
        result.solution.allFinite() &&
        result.residual <= tolerance;
    return result;
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
        case PressureSolver2D::Jacobi: {
            return SolveStationary(
                StationaryPressureSolver2D::Jacobi,
                matrix,
                rhs,
                max_iterations,
                tolerance);
        }
        case PressureSolver2D::GaussSeidel: {
            return SolveStationary(
                StationaryPressureSolver2D::GaussSeidel,
                matrix,
                rhs,
                max_iterations,
                tolerance);
        }
        case PressureSolver2D::SuccessiveOverRelaxation: {
            return SolveStationary(
                StationaryPressureSolver2D::SuccessiveOverRelaxation,
                matrix,
                rhs,
                max_iterations,
                tolerance);
        }
        case PressureSolver2D::ConjugateGradient: {
            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<Real>,
                Eigen::Lower | Eigen::Upper,
                Eigen::DiagonalPreconditioner<Real>> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::ConjugateGradientNoPreconditioner: {
            Eigen::ConjugateGradient<
                Eigen::SparseMatrix<Real>,
                Eigen::Lower | Eigen::Upper,
                Eigen::IdentityPreconditioner> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::BiCGSTAB: {
            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::LeastSquaresConjugateGradient: {
            Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
        }
        case PressureSolver2D::SimplicialLLT: {
            Eigen::SimplicialLLT<Eigen::SparseMatrix<Real>> eigen_solver;
            return SolveDirect(eigen_solver, matrix, rhs);
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
        Eigen::Lower | Eigen::Upper,
        Eigen::DiagonalPreconditioner<Real>> eigen_solver;
    return SolveIterative(eigen_solver, matrix, rhs, max_iterations, tolerance);
}

}  // namespace skyspaces
