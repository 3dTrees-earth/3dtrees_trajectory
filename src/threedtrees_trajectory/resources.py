"""Process resource accounting and the CPU budget of a Galaxy or local job."""
from __future__ import annotations

import os
import resource
import sys
from pathlib import Path
from time import perf_counter, process_time

CGROUP = Path("/sys/fs/cgroup")


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="ascii").strip()
    except OSError:
        return ""


def _cgroup() -> Path | None:
    for line in _read(Path("/proc/self/cgroup")).splitlines():
        if line.startswith("0::/"):
            # A private cgroup namespace (as in Docker) exposes its root as '/'.
            group = CGROUP / line[4:]
            if group.is_dir() and group.is_relative_to(CGROUP) and ".." not in group.parts:
                return group
    return None


def _limit(value: str) -> int | None:
    return int(value) if value.isascii() and value.isdecimal() else None


def decoding_threads(requested: int | None) -> int:
    """Limit work to Galaxy slots, CPU affinity, and any ancestor cgroup v2 quota."""
    configured = _limit(os.environ.get("GALAXY_SLOTS", "")) or None
    budgets = [requested or configured or 8]
    if configured is not None:
        budgets.append(configured)
    if hasattr(os, "sched_getaffinity"):
        budgets.append(len(os.sched_getaffinity(0)))
    else:
        budgets.append(os.cpu_count() or 1)
    group = _cgroup()
    if group is not None:
        for path in (group, *group.parents):
            if not path.is_relative_to(CGROUP):
                break
            quota, _, period = _read(path / "cpu.max").partition(" ")
            if _limit(quota) is not None and _limit(period):
                budgets.append(max(1, int(quota) // int(period)))
    return min(budgets)


class RunResources:
    def __init__(self) -> None:
        self.start_wall = perf_counter()
        self.start_cpu = process_time()
        self.group = _cgroup()
        self.start_group_cpu = self._group_cpu()

    def _group_cpu(self) -> int | None:
        if self.group is not None:
            for line in _read(self.group / "cpu.stat").splitlines():
                key, _, value = line.partition(" ")
                if key == "usage_usec":
                    return _limit(value)
        return None

    def report(self, threads: int | None, points_read: int) -> dict:
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        data = {
            "wall_seconds": round(perf_counter() - self.start_wall, 6),
            "cpu_seconds": round(process_time() - self.start_cpu, 6),
            "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
            "points_read": points_read,
            "decoding_threads": threads,
        }
        if self.group is not None:
            end_cpu = self._group_cpu()
            cpu_seconds = None
            if end_cpu is not None and self.start_group_cpu is not None:
                cpu_seconds = (end_cpu - self.start_group_cpu) / 1_000_000
            data["cgroup_v2"] = {
                "memory_current_bytes": _limit(_read(self.group / "memory.current")),
                "memory_peak_bytes": _limit(_read(self.group / "memory.peak")),
                "memory_limit_bytes": _limit(_read(self.group / "memory.max")),
                "cpu_seconds": cpu_seconds,
            }
        return data
