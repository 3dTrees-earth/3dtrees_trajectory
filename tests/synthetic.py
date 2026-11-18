"""Synthetic lidar clouds with a known sensor: rays from the sensor hit a 'canopy' surface (first return) and the
ground (last return)."""
from pathlib import Path

import laspy
import numpy as np


def _write(path: Path, xyz, time, return_number, number_of_returns, scan_angle, source=None, *,
           point_format=6, version="1.4"):
    header = laspy.LasHeader(point_format=point_format, version=version)
    header.offsets = np.floor(xyz.min(axis=0))
    header.scales = [0.001, 0.001, 0.001]
    las = laspy.LasData(header)
    las.x, las.y, las.z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    if "gps_time" in header.point_format.dimension_names:
        las.gps_time = time
    las.return_number = return_number
    las.number_of_returns = number_of_returns
    if point_format >= 6:
        las.scan_angle = np.round(scan_angle / 0.006).astype(np.int16)
    else:
        las.scan_angle_rank = np.round(scan_angle).astype(np.int8)
    if source is not None:
        las.point_source_id = source
    las.write(path)


def _hits(sensor, direction, rng, canopy_share):
    """First hit on a canopy layer (for canopy_share of pulses) at 10-20 m, last hit on the ground z = 0."""
    to_ground = -sensor[:, 2] / direction[:, 2]
    ground = sensor + to_ground[:, None] * direction
    canopy_z = rng.uniform(10, 20, len(sensor))
    to_canopy = (canopy_z - sensor[:, 2]) / direction[:, 2]
    canopy = sensor + to_canopy[:, None] * direction
    has_canopy = rng.random(len(sensor)) < canopy_share
    return canopy, ground, has_canopy


def _records(time, first, last, has_two, angle, single_return):
    if single_return:
        start = np.where(has_two[:, None], first, last)
        n = len(time)
        return start, time, np.ones(n, np.uint8), np.ones(n, np.uint8), angle
    xyz, t, rn, nr, ang = [], [], [], [], []
    for i in range(len(time)):  # small synthetic sets only
        if has_two[i]:
            xyz.extend((first[i], last[i]))
            t.extend((time[i], time[i]))
            rn.extend((1, 2))
            nr.extend((2, 2))
            ang.extend((angle[i], angle[i]))
        else:
            xyz.append(last[i])
            t.append(time[i])
            rn.append(1)
            nr.append(1)
            ang.append(angle[i])
    return np.array(xyz), np.array(t), np.array(rn, np.uint8), np.array(nr, np.uint8), np.array(ang)


def flight(path: Path, single_return=False, seconds=20.0, rate=4000, seed=0, *,
           point_format=6, version="1.4"):
    """Line scanner flying north at 5 m/s, 60 m above ground, scanning +-40 degrees across track."""
    rng = np.random.default_rng(seed)
    time = 300000.0 + np.arange(0, seconds, 1 / rate)
    sensor = np.c_[np.full(len(time), 1000.0), 2000.0 + 5.0 * (time - time[0]), np.full(len(time), 60.0)]
    angle = rng.uniform(-40, 40, len(time))
    a = np.radians(angle)
    direction = np.c_[np.sin(a), np.zeros(len(a)), -np.cos(a)]
    first, last, two = _hits(sensor, direction, rng, 0.6)
    _write(path, *_records(time, first, last, two, angle, single_return),
           point_format=point_format, version=version)
    return time, sensor


def scan(path: Path, position, single_return=False, rays=40000, seed=0, t0=100000.0, source=None):
    """Static scanner at position; rays fan out below the horizon onto canopy and ground."""
    rng = np.random.default_rng(seed)
    azimuth = rng.uniform(0, 2 * np.pi, rays)
    zenith = np.radians(rng.uniform(92, 175, rays))  # downward rays only (out to about 40 m), so the ground catches each one
    direction = np.c_[np.sin(zenith) * np.cos(azimuth), np.sin(zenith) * np.sin(azimuth), np.cos(zenith)]
    sensor = np.tile(np.asarray(position, float), (rays, 1))
    # 'canopy' here: low shrubs 0.2-1 m above ground between sensor and ground
    to_ground = -sensor[:, 2] / direction[:, 2]
    ground = sensor + to_ground[:, None] * direction
    frac = rng.uniform(0.3, 0.8, rays)
    shrub = sensor + (frac * to_ground)[:, None] * direction
    two = rng.random(rays) < 0.5
    time = t0 + np.arange(rays) * 1e-4
    xyz, t, rn, nr, ang = _records(time, shrub, ground, two, np.zeros(rays), single_return)
    _write(path, xyz, t, rn, nr, ang, None if source is None else np.full(len(t), source, np.uint16))
