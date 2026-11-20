// Voxel stepping and echo classification follow OccPy 0.3.0 (MIT), Raytracer.cpp.
// OccPy Raytracer.cpp copyright 2023 Daniel Kükenbrink.
// See https://github.com/dkueken/OccPy. Count storage is limited to relevant
// columns, but entry, voxel faces and step ties use the whole acquisition grid.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {
using XYZ = std::array<double, 3>;
using Cell = std::array<int32_t, 3>;

struct Grid {
    const double* lower;
    const double* upper;
    const int32_t* shape;

    Cell cell(const XYZ& point) const {
        Cell result;
        for (int axis = 0; axis < 3; ++axis) {
            const auto size = upper[axis] - lower[axis];
            result[axis] = static_cast<int32_t>(std::floor((point[axis] - lower[axis]) / size * shape[axis]));
        }
        return result;
    }

    bool contains(const Cell& point) const {
        for (int axis = 0; axis < 3; ++axis) {
            if (point[axis] < 0 || point[axis] >= shape[axis]) {
                return false;
            }
        }
        return true;
    }

    bool entry(const double* sensor, const XYZ& direction, double& enter) const {
        double leave = std::numeric_limits<double>::infinity();
        enter = -leave;
        for (int axis = 0; axis < 3; ++axis) {
            const auto a = (lower[axis] - sensor[axis]) / direction[axis];
            const auto b = (upper[axis] - sensor[axis]) / direction[axis];
            const auto near = std::min(a, b);
            const auto far = std::max(a, b);
            if (enter > far || near > leave) {
                return false;
            }
            enter = std::max(enter, near);
            leave = std::min(leave, far);
        }
        enter = std::max(0.0, enter);
        return true;
    }
};

struct Counts {
    const int32_t* rectangle;  // x start, y start, x stop, y stop
    const int32_t* first_layer;
    const int64_t* offsets;
    uint32_t* values;

    uint32_t* at(const Cell& point) const {
        const auto x = point[0], y = point[1];
        if (x < rectangle[0] || y < rectangle[1] || x >= rectangle[2] || y >= rectangle[3]) {
            return nullptr;
        }
        const int64_t column = int64_t(x - rectangle[0]) * (rectangle[3] - rectangle[1]) + y - rectangle[1];
        const int64_t depth = point[2] - first_layer[column];
        const int64_t start = offsets[column];
        if (depth < 0 || depth >= offsets[column + 1] - start) {
            return nullptr;
        }
        return values + 3 * (start + depth);
    }
};

double distance_squared(const XYZ& point, const double* origin) {
    double result = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        const auto delta = point[axis] - origin[axis];
        result += delta * delta;
    }
    return result;
}

void trace_pulse(const Grid& grid, const Counts& counts, const double* coords, const double* sensor,
                 const uint8_t* returns, int n, int64_t* stats) {
    std::array<XYZ, 15> echoes;
    for (int i = 0; i < n; ++i) {
        auto& echo = echoes[returns ? returns[i] - 1 : 0];
        std::copy_n(coords + 3 * i, 3, echo.begin());
    }
    XYZ direction;
    for (int axis = 0; axis < 3; ++axis) {
        direction[axis] = echoes[n - 1][axis] - sensor[axis];
        if (!std::isfinite(direction[axis])) {
            return;
        }
        if (direction[axis] == 0.0) {
            direction[axis] = std::numeric_limits<double>::min();
        }
    }
    double enter;
    if (!grid.entry(sensor, direction, enter) || !std::isfinite(enter)) {
        return;
    }
    ++stats[0];
    XYZ start;
    for (int axis = 0; axis < 3; ++axis) {
        start[axis] = sensor[axis] + enter * direction[axis];
    }

    std::array<Cell, 15> targets;
    int outside = 0;
    int hit = 0;
    const auto entry_distance = std::sqrt(distance_squared(start, sensor));
    for (int i = 0; i < n; ++i) {
        targets[i] = grid.cell(echoes[i]);
        outside += !grid.contains(targets[i]);
        hit += std::sqrt(distance_squared(echoes[i], sensor)) < entry_distance;
    }
    Cell position = grid.cell(start);
    Cell step;
    XYZ face, delta;
    for (int axis = 0; axis < 3; ++axis) {
        if (position[axis] == grid.shape[axis]) {
            --position[axis];
        }
        const bool forward = direction[axis] >= 0.0;
        step[axis] = forward ? 1 : -1;
        const auto size = grid.upper[axis] - grid.lower[axis];
        const auto fraction = double(position[axis] + int(forward)) / grid.shape[axis];
        const auto boundary = grid.lower[axis] + fraction * size;
        face[axis] = enter + (boundary - start[axis]) / direction[axis];
        delta[axis] = (size / grid.shape[axis]) / std::abs(direction[axis]);
    }

    while (grid.contains(position)) {
        auto* cell = counts.at(position);
        if (hit < n && targets[hit] == position) {
            do {
                if (cell) {
                    ++cell[0];
                }
                ++hit;
            } while (hit < n && targets[hit] == position);
        } else if (cell) {
            ++cell[hit < n ? 1 : 2];
        }
        const auto next = std::min({face[0], face[1], face[2]});
        for (int axis = 0; axis < 3; ++axis) {
            if (face[axis] == next) {
                position[axis] += step[axis];
                face[axis] += delta[axis];
            }
        }
    }
    if (returns && hit < n - outside) {
        stats[1] += n - hit;
    }
}
}  // namespace

extern "C" void trace_rays(const double* lower, const double* upper, const int32_t* shape,
                           const int32_t* rectangle, const int32_t* first_layer, const int64_t* offsets,
                           uint32_t* values, const double* coords, const double* sensors,
                           const uint8_t* returns, const int32_t* pulses, int32_t n, int64_t* stats) {
    const Grid grid{lower, upper, shape};
    const Counts counts{rectangle, first_layer, offsets, values};
    for (int32_t i = 0; i < n; ++i) {
        const int32_t first = pulses ? pulses[i] : i;
        const int32_t size = pulses ? pulses[i + 1] - first : 1;
        trace_pulse(grid, counts, coords + 3LL * first, sensors + 3LL * i,
                    returns ? returns + first : nullptr, size, stats);
    }
}

extern "C" void sum_counts(const int32_t* rectangle, const int32_t* first_layer, const int64_t* offsets,
                            uint32_t* source, const int32_t* start, const int32_t* stop, uint64_t* dest) {
    const int64_t ny = stop[1] - start[1], nz = stop[2] - start[2];
    for (int32_t x = std::max(start[0], rectangle[0]); x < std::min(stop[0], rectangle[2]); ++x) {
        for (int32_t y = std::max(start[1], rectangle[1]); y < std::min(stop[1], rectangle[3]); ++y) {
            const int64_t col = int64_t(x - rectangle[0]) * (rectangle[3] - rectangle[1]) + y - rectangle[1];
            const int64_t base = offsets[col];
            const auto first = first_layer[col];
            const int64_t lower = std::max(start[2], first);
            const int64_t upper = std::min(int64_t(stop[2]), first + offsets[col + 1] - base);
            if (lower >= upper) {
                continue;
            }
            auto* to = dest + 3 * (((x - start[0]) * ny + (y - start[1])) * nz + lower - start[2]);
            const auto* from = source + 3 * (base + lower - first);
            for (int64_t i = 0; i < 3 * (upper - lower); ++i) {
                to[i] += from[i];
            }
        }
    }
}
