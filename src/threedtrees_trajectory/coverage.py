"""Reasons why a moving acquisition cannot produce a trajectory or occlusion."""
from __future__ import annotations


class UnusableTime(ValueError):
    exit_code = 3


class NoTrajectory(ValueError):
    exit_code = 4
