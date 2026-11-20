"""Whole-grid ray stepping with compact counts for the columns used in the products."""
from __future__ import annotations

import ctypes
import os
import sysconfig
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np


@lru_cache
def library():
    directory = Path(__file__).parent
    names = (f"_ray_walk{sysconfig.get_config_var('EXT_SUFFIX')}", "_ray_walk.so")
    candidates = [directory / name for name in names]
    if location := os.environ.get("TRAJECTORY_RAY_WALK"):
        candidates.insert(0, Path(location))
    for candidate in candidates:
        if candidate.is_file():
            lib = ctypes.CDLL(str(candidate))
            lib.trace_rays.argtypes = [ctypes.c_void_p] * 11 + [ctypes.c_int32, ctypes.c_void_p]
            lib.trace_rays.restype = None
            lib.sum_counts.argtypes = [ctypes.c_void_p] * 7
            lib.sum_counts.restype = None
            return lib
    raise ImportError("The trajectory ray walker is missing. Install the package or run the tool image.")


@dataclass(frozen=True)
class Columns:
    rectangle: np.ndarray
    first: np.ndarray
    offsets: np.ndarray

    @property
    def cells(self):
        return int(self.offsets[-1])

    @classmethod
    def from_bounds(cls, bounds, first, stop):
        xs, ys, zs = bounds
        rectangle = np.array([xs.start, ys.start, xs.stop, ys.stop], np.int32)
        first_layer = np.ascontiguousarray(np.maximum(first[xs, ys], zs.start), dtype=np.int32).ravel()
        final_layer = np.minimum(stop[xs, ys], zs.stop).ravel()
        lengths = np.maximum(0, final_layer - first_layer)
        offsets = np.empty(len(lengths) + 1, np.int64)
        offsets[0] = 0
        np.cumsum(lengths, out=offsets[1:])
        return cls(rectangle, first_layer, offsets)


def layers(grid, top, ground):
    """Keep the band from the ground to the highest raw echo plus its possible face nudge."""
    nz = grid.shape[2]
    first = np.floor((ground - grid.origin[2]) / grid.voxel) - 1
    first = np.clip(first, 0, nz).astype(np.int32)
    stop = np.where(top >= 0, np.minimum(top + 2, nz), 0).astype(np.int32)
    return first, stop


class Walker:
    def __init__(self, grid, columns, counts):
        self.grid_args = (
            np.ascontiguousarray(grid.origin, dtype=np.float64),
            np.ascontiguousarray(grid.origin + np.array(grid.shape) * grid.voxel, dtype=np.float64),
            np.asarray(grid.shape, dtype=np.int32),
        )
        self.count_args = columns.rectangle, columns.first, columns.offsets, counts
        self.stats = np.zeros(2, np.int64)
        self.trace = library().trace_rays

    def add(self, coordinates, sensors, returns=None, pulses=None):
        coords, xyz = np.column_stack(coordinates), np.column_stack(tuple(sensors))
        self.trace(*(array.ctypes.data for array in (*self.grid_args, *self.count_args, coords, xyz)),
                   returns.ctypes.data if returns is not None else None,
                   pulses.ctypes.data if pulses is not None else None, len(xyz), self.stats.ctypes.data)


def sum_block(columns, values, start, stop, counts):
    library().sum_counts(*(array.ctypes.data for array in
                           (columns.rectangle, columns.first, columns.offsets, values, start, stop, counts)))
