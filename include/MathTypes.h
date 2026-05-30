#pragma once

#include <Eigen/Core>

namespace skyspaces {

using Real = double;
using Vector2R = Eigen::Matrix<Real, 2, 1>;
using Vector2I = Eigen::Vector2i;

// Scalar grid storage. Rows are x and columns are y so grid(x, y) maps
// directly to data(x, y). Row-major layout keeps y-contiguous inner loops fast.
using ScalarArray2D = Eigen::Array<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

}  // namespace skyspaces
