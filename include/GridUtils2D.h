#pragma once

#include "MathTypes.h"

#include <algorithm>
#include <cassert>

namespace skyspaces {

inline int FlattenCellIndex2D(int x, int y, int height) {
    assert(height > 0);
    return x * height + y;
}

inline int ClampGridIndex2D(int value, int size) {
    assert(size > 0);
    return std::clamp(value, 0, size - 1);
}

inline Real CellCenterGridCoordinate2D(Real world_coordinate, Real inverse_cell_size) {
    return world_coordinate * inverse_cell_size - 0.5;
}

inline Real FaceGridCoordinate2D(Real world_coordinate, Real inverse_cell_size) {
    return world_coordinate * inverse_cell_size;
}

}  // namespace skyspaces
