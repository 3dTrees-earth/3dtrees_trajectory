"""Render the stage and decision figures from the checked-in measurements."""
import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output", type=Path)
parser.add_argument("--stages", type=Path, default=Path(__file__).with_name("stages.csv"))
parser.add_argument("--scaling", type=Path, help="Trajectory release benchmark CSV")
args = parser.parse_args()
OUT = args.output
OUT.mkdir(parents=True, exist_ok=True)
with args.stages.open() as stream:
    table = list(csv.DictReader(stream))

stages = [
    ("Read + path", ("read/time-sort", "platform check", "path estimation")),
    ("Extent + routing", ("grid extent/crop pass", "route tile spools", "finalize extent/crop", "plan GPS-time windows")),
    ("Traverse", ("traversal",)),
    ("Merge + write", ("classify + ground", "classify voxels + write products", "COPC write", "COG write")),
]
colors = ["#4878a8", "#d9822b", "#34846d", "#8256a2"]
rows = {(r["run"], r["mode"], r["stage"]): r for r in table}


def plot_stages(modes, title, filename, *, minutes=False):
    runs = list(dict.fromkeys((r["run"], r["mode"]) for r in table if r["mode"] in modes))
    if not runs:
        return
    labels = [f"{name} · {mode}" for name, mode in runs]
    fig, axes = plt.subplots(1, 3, figsize=(14, max(3.5, 0.65 * len(runs) + 1.5)),
                             gridspec_kw={"width_ratios": [1.8, 1, 1]})
    positions = range(len(runs))
    divisor = 60 if minutes else 1
    totals = []
    for position, run in enumerate(runs):
        offset = 0.0
        for (name, members), color in zip(stages, colors):
            seconds = sum(float(rows[(*run, member)]["wall_s"]) for member in members if (*run, member) in rows)
            duration = seconds / divisor
            axes[0].barh(position, duration, left=offset, color=color, label=name if position == 0 else None)
            offset += duration
        totals.append(offset)
        read = rows.get((*run, "read/time-sort"), {})
        trace = rows.get((*run, "traversal"), {})
        axes[1].barh(position - 0.15, float(read.get("cores", 0)), height=0.3, color=colors[0],
                     label="Read" if position == 0 else None)
        axes[1].barh(position + 0.15, float(trace.get("cores", 0)), height=0.3, color=colors[2],
                     label="Traversal" if position == 0 else None)
        peak = max(float(r["peak_ram_gib"]) for r in table if (r["run"], r["mode"]) == run)
        axes[2].barh(position, peak, color="#748399")

    width = max(totals)
    legend_y = -0.25 if minutes else -0.16
    for position, duration in enumerate(totals):
        axes[0].text(duration + width * 0.02, position, f"{duration:.0f}", va="center", fontsize=8)
    axes[0].set(title="Wall time by stage", xlabel="Minutes" if minutes else "Seconds",
                yticks=list(positions), yticklabels=labels, xlim=(0, width * 1.16))
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, legend_y), ncol=2, frameon=False)
    axes[1].set(title="Read and traversal CPU", xlabel="CPU seconds / wall second", xlim=(0, 10.5), yticks=[])
    axes[1].axvline(10, color="#555", linestyle="--", linewidth=0.8)
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, legend_y), ncol=2, frameon=False)
    axes[2].set(title="Peak stage RSS", xlabel="GiB, process or sampled worker sum", yticks=[])
    for ax in axes:
        ax.set_ylim(len(runs) - 0.5, -0.5)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(title, y=0.98)
    fig.tight_layout()
    fig.savefig(OUT / filename, dpi=160)
    plt.close(fig)


plot_stages(("Round 9 control", "Round 10"), "Trajectory controls · 10 CPU slots", "stages.png")
plot_stages(("Round 8 50 GiB", "Round 10 50 GiB"), "Dataset 2585 full replays · 10 CPU slots and 50 GiB",
            "stages-2585.png", minutes=True)

def plot_scaling(path):
    with path.open(newline="") as stream:
        measurements = list(csv.DictReader(stream))

    slots = (2, 4, 6, 8, 10)
    datasets = (("1807", "1807 ULS"), ("2014", "2014 ULS"), ("1844", "1844 ULS"), ("1788", "1788 MLS"))
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), layout="constrained")
    for (dataset, label), color, marker in zip(datasets, colors, ("o", "s", "^", "D")):
        runs = sorted((row for row in measurements if row["dataset"] == dataset),
                      key=lambda row: int(row["threads"]))
        assert tuple(int(row["threads"]) for row in runs) == slots, dataset
        speed = [float(row["speedup_vs_10"]) for row in runs]
        wall_ratio = [1 / value for value in speed]
        allocation = [n * value / 10 for n, value in zip(slots, wall_ratio)]
        memory = [float(row["peak_anon_gib"]) for row in runs]
        for ax, values in zip(axes, (wall_ratio, allocation, memory)):
            ax.plot(slots, values, color=color, marker=marker, label=label)

    for ax, title, units in zip(
        axes,
        ("Elapsed time", "CPU allocation", "Peak anonymous RSS"),
        ("× 10-slot elapsed time", "× 10-slot core-seconds", "GiB"),
    ):
        ax.set(title=title, xlabel="CPU slots", ylabel=units, xticks=slots)
        ax.grid(axis="y", alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[:2]:
        ax.axhline(1, color="#555", linestyle="--", linewidth=0.8, zorder=0)
    axes[0].legend(loc="upper right", frameon=False, fontsize=8)
    fig.suptitle("Trajectory 1.0.0 CPU scaling · 50 GiB per job", fontsize=14)
    fig.savefig(OUT / "cpu-scaling.png", dpi=160)
    plt.close(fig)


if args.scaling is not None:
    plot_scaling(args.scaling)


# ---------- flow.png: decision flow ----------
fig, ax = plt.subplots(figsize=(13, 15))
ax.set_xlim(0, 13)
ax.set_ylim(0, 30)
ax.axis("off")


def box(x, y, text, color="#e8f0fe", w=4.6, h=1.5):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.08", fc=color, ec="#333", lw=1))
    ax.text(x, y, text, ha="center", va="center", fontsize=9.5, wrap=True)


def stop(x, y, text):
    box(x, y, text, color="#fde2e1", w=3.6, h=1.3)


def ok(x, y, text):
    box(x, y, text, color="#e3f5e1", w=4.2, h=1.4)


def arrow(x1, y1, x2, y2, label=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="->", lw=1.1, color="#333"))
    if label:
        ax.text((x1 + x2) / 2 + 0.1, (y1 + y2) / 2 + 0.1, label, fontsize=8.5, color="#333")


box(5, 29, "Input: one collection of standardized LAS/LAZ\n(+ optional delivered trajectory CSV and terrain raster)")
box(5, 26.6, "Argument + path checks\n(outputs ≠ inputs, voxel 0.05–10 m, ≥ 1 file)")
stop(10.6, 26.6, "exit 2: invalid\narguments / unreadable input")
box(5, 24.2, "Read scanner-order chunks and index bounds\nGPS time present, finite, varying?")
stop(10.6, 24.2, "exit 3: no / constant /\nnon-finite GPS time")
box(5, 21.8, "Rebuild pulses (return pairing per file)\nusable pulses?")
stop(10.6, 21.8, "exit 4: no usable\npulses")
box(5, 19.4, "Platform check: static scanner (TLS)?\n(lines meet ≤ 5 cm, or ground density ≥ 15×)")
stop(10.6, 19.4, "exit 2: TLS refused\n(v1.1 on branch)")
box(5, 16.8, "Path method per file\n1 pulse lines (multi-return)  2 scan angle (single)\n3 lines + scan-angle gap fill", w=5.4, h=1.9)
stop(10.6, 16.8, "exit 4: single returns\nwithout scan angle,\nno delivered path")
box(5, 13.9, "Quality: median line miss ≤ 0.1 m good · ≤ 0.5 m rough\n> 0.5 m unreliable → retry with scan angle", w=5.6, h=1.6)
ok(2.3, 11.4, "Trajectory CSV\nJSON report with --debug")
box(7.4, 11.4, "Occlusion requested?\npath = delivered (always, if given;\nflag if median miss > 0.1 m) else estimated", w=5.2, h=1.8)
box(7.4, 8.6, "Merge first-read column bounds in path coverage\nSample terrain raster with lowest-hit fallback", w=5.6, h=1.6)
box(7.4, 6.0, "Stream disjoint GPS windows in parallel\nAdjust admission using observed worker peak RSS.", w=5.6, h=1.6)
box(7.4, 3.5, "Classify blocks in parallel, then write in grid order", w=5.8, h=1.2)
ok(7.4, 1.2, "Occlusion COPC (classes 1, 3 and 4)\n+ COG (one band per voxel layer, 0.5 m default)")
box(11.65, 8.6, "Occlusion error\nstderr, empty products,\nexit 2 if scratch is full", color="#fde2e1", w=2.4, h=1.9)

for y1, y2 in ((28.2, 27.4), (25.8, 25.0), (23.4, 22.6), (21.0, 20.2), (18.6, 17.8), (15.8, 14.7)):
    arrow(5, y1, 5, y2)
for y in (26.6, 24.2, 21.8, 19.4, 16.8):
    arrow(7.3 if y != 16.8 else 7.7, y, 8.8, y, "no" if y in (24.2, 21.8) else ("yes" if y == 19.4 else None))
arrow(4.2, 13.1, 2.6, 12.1)
arrow(5.8, 13.1, 7.0, 12.3)
arrow(7.4, 10.5, 7.4, 9.4)
arrow(7.4, 7.8, 7.4, 6.8)
arrow(7.4, 5.2, 7.4, 4.1)
arrow(7.4, 2.9, 7.4, 2.0)
arrow(10.25, 8.6, 10.4, 8.6)
ax.text(0.2, 0.2, "Exit codes: 0 done · 2 invalid input / TLS · 3 GPS time unusable · 4 no trajectory possible. "
        "With --debug, refusals write a JSON report.", fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "flow.png", dpi=150)
plt.close(fig)
