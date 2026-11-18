"""Process resource accounting and the CPU budget of a Galaxy or local job."""
from __future__ import annotations

import os
import resource
import sys
from contextlib import contextmanager
from pathlib import Path
from threading import Event, Thread
from time import perf_counter

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


def _cpu_seconds() -> float:
    return sum(usage.ru_utime + usage.ru_stime for usage in
               (resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)))


def _rss(kind=resource.RUSAGE_SELF) -> int:
    value = resource.getrusage(kind).ru_maxrss
    return value if sys.platform == "darwin" else value * 1024


def _proc_rss(pid: int | str, key: str) -> int:
    for line in _read(Path(f"/proc/{pid}/status")).splitlines():
        if line.startswith(f"{key}:"):
            return int(line.split()[1]) * 1024
    return 0


class RunResources:
    def __init__(self) -> None:
        self.start_wall = perf_counter()
        self.start_cpu = _cpu_seconds()
        self.group = _cgroup()
        self.start_group_cpu = self._group_cpu()
        self.peak = 0
        self.stages = []

    @contextmanager
    def stage(self, name: str, *, children: bool = False):
        start_wall, start_cpu = perf_counter(), _cpu_seconds()
        self.peak = max(self.peak, _rss())
        try:
            Path("/proc/self/clear_refs").write_text("5", encoding="ascii")
            reset = True
        except OSError:
            reset = False
        child_peak_start = _rss(resource.RUSAGE_CHILDREN)
        child_peak = 0
        stop = Event()

        def monitor_children():
            nonlocal child_peak
            tasks = Path("/proc/self/task")
            while not stop.is_set():
                if tasks.is_dir():
                    for task in tasks.iterdir():
                        for pid in _read(task / "children").split():
                            child_peak = max(child_peak, _proc_rss(pid, "VmHWM"))
                stop.wait(0.1)

        monitor = Thread(target=monitor_children, daemon=True) if children else None
        if monitor is not None:
            monitor.start()
        stage = {"name": name}
        try:
            yield stage
        finally:
            if monitor is not None:
                stop.set()
                monitor.join()
            peak = (_proc_rss("self", "VmHWM") if reset else 0) or _rss()
            reaped_peak = _rss(resource.RUSAGE_CHILDREN)
            if children and reaped_peak > child_peak_start:
                child_peak = max(child_peak, reaped_peak)
            self.peak = max(self.peak, peak, child_peak)
            stage.update(wall_seconds=round(perf_counter() - start_wall, 6),
                         cpu_seconds=round(_cpu_seconds() - start_cpu, 6),
                         peak_rss_bytes=max(peak, child_peak), rss_scope="stage" if reset else "process-lifetime")
            if children:
                stage["child_peak_rss_bytes"] = child_peak
            self.stages.append(stage)

    def _group_cpu(self) -> int | None:
        if self.group is not None:
            for line in _read(self.group / "cpu.stat").splitlines():
                key, _, value = line.partition(" ")
                if key == "usage_usec":
                    return _limit(value)
        return None

    def report(self, threads: int | None, points_read: int) -> dict:
        data = {
            "wall_seconds": round(perf_counter() - self.start_wall, 6),
            "cpu_seconds": round(_cpu_seconds() - self.start_cpu, 6),
            "peak_rss_bytes": max(self.peak, _rss()),
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
