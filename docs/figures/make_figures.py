"""Build the stage table and figures for docs/REPORT.md from measured run records (rounds 2-3).

Input: /mnt/ssds/kg281/runs/2026-11-20_trajectory_codex-round3/stage-resources.tsv (stage totals = rows without tile).
Output: stages.csv, stages.png, flow.png in the folder given as argv[1].
"""
import csv
import sys
from collections import OrderedDict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
TSV = Path("/mnt/ssds/kg281/runs/2026-11-20_trajectory_codex-round3/stage-resources.tsv")

RUNS = OrderedDict([  # label -> (acquisition key, points, input note)
    ("1807 ULS\n53 M", ("1807 LAZ", 52_718_335, "LAZ, round 2")),
    ("2014 ULS\nmulti", ("2014 multi LAZ", None, "LAZ, round 2")),
    ("1788 MLS\n25 M wide", ("mls-occlusion/1788", 24_653_752, "standardized LAZ, round 3")),
    ("1844 ULS\n678 M", ("1844 LAZ", 678_104_137, "LAZ, round 2")),
    ("2585 ULS\n7.0 B", ("2585 LAZ", 7_044_924_732, "LAZ, round 2")),
])
STAGES = ["read/time-sort", "platform check", "path estimation", "grid extent/crop pass", "route tile spools",
          "traversal", "classify + ground", "COPC write", "COG write"]
SHORT = {"read/time-sort": "read + time sort", "grid extent/crop pass": "extent / crop pass",
         "route tile spools": "route to tiles", "classify + ground": "classify + ground"}

rows = list(csv.DictReader(TSV.open(), delimiter="\t"))


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


table = []
for label, (key, points, note) in RUNS.items():
    for stage in STAGES:
        hits = [r for r in rows if r["acquisition"] == key and r["name"] == stage and not r["tile"]]
        if not hits:
            continue
        r = hits[-1]
        wall, cpu = num(r["wall_seconds"]), num(r["cpu_seconds"])
        rss = max(v for v in (num(r["peak_rss_bytes"]), num(r["sampled_aggregate_rss_bytes"]), 0.0) if v is not None)
        scratch = num(r["sampled_peak_scratch_bytes"])
        table.append(dict(run=label.replace("\n", " "), input=note, stage=stage, wall_s=round(wall or 0, 1),
                          cores=round(cpu / wall, 1) if wall and cpu else None, peak_ram_gib=round(rss / 2**30, 1),
                          scratch_gib=None if scratch is None else round(scratch / 2**30, 1)))

with (OUT / "stages.csv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(table[0]))
    w.writeheader()
    w.writerows(table)

# ---------- stages.png: wall share, cores used, peak RAM ----------
colors = plt.cm.tab10.colors
fig, axes = plt.subplots(1, 3, figsize=(17, 6.2), gridspec_kw=dict(width_ratios=[1.2, 1, 1]))
labels = list(RUNS)
for i, stage in enumerate(STAGES):
    left = []
    for label in labels:
        name = label.replace("\n", " ")
        prior = sum(t["wall_s"] for t in table if t["run"] == name and STAGES.index(t["stage"]) < i)
        left.append(prior)
    walls = [next((t["wall_s"] for t in table if t["run"] == l.replace("\n", " ") and t["stage"] == stage), 0)
             for l in labels]
    totals = [sum(t["wall_s"] for t in table if t["run"] == l.replace("\n", " ")) for l in labels]
    share = [w / tot * 100 if tot else 0 for w, tot in zip(walls, totals)]
    start = [p / tot * 100 if tot else 0 for p, tot in zip(left, totals)]
    axes[0].barh(range(len(labels)), share, left=start, color=colors[i % 10], label=SHORT.get(stage, stage))
for j, label in enumerate(labels):
    total = sum(t["wall_s"] for t in table if t["run"] == label.replace("\n", " "))
    axes[0].text(101, j, f"{total / 3600:.2f} h" if total > 3600 else f"{total / 60:.0f} min", va="center", fontsize=9)
axes[0].set_yticks(range(len(labels)), labels)
axes[0].invert_yaxis()
axes[0].set_xlim(0, 115)
axes[0].set_xlabel("share of wall time (%)  ·  total on the right")
axes[0].set_title("Where the time goes")
axes[0].legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=3, frameon=False)

key_stages = ["read/time-sort", "path estimation", "grid extent/crop pass", "route tile spools", "traversal",
              "COPC write", "COG write"]
width = 0.15
for k, label in enumerate(labels):
    name = label.replace("\n", " ")
    cores = [next((t["cores"] or 0 for t in table if t["run"] == name and t["stage"] == s), 0) for s in key_stages]
    ram = [next((t["peak_ram_gib"] for t in table if t["run"] == name and t["stage"] == s), 0) for s in key_stages]
    x = [i + (k - 2) * width for i in range(len(key_stages))]
    axes[1].bar(x, cores, width, label=name, color=colors[k])
    axes[2].bar(x, ram, width, label=name, color=colors[k])
for ax, title, ylabel in ((axes[1], "CPU cores used (of 10)", "cores = CPU s / wall s"),
                          (axes[2], "Peak RAM per stage", "GiB (largest process or sampled sum of workers)")):
    ax.set_xticks(range(len(key_stages)), [SHORT.get(s, s) for s in key_stages], rotation=35, ha="right", fontsize=9)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
axes[1].axhline(10, color="grey", lw=0.8, ls="--")
axes[1].legend(fontsize=8)
fig.suptitle("3dtrees_trajectory: measured stage resources (LAZ input, 10-core cap; 2585 COG write = old pixel-interleaved "
             "intermediate, now 51 s)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "stages.png", dpi=150)

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


box(5, 29, "Input: one collection of standardized LAS/LAZ\n(+ optional delivered trajectory CSV, optional DTM)")
box(5, 26.6, "Argument + path checks\n(outputs ≠ inputs, voxel 0.05–10 m, ≥ 1 file)")
stop(10.6, 26.6, "exit 2: invalid\narguments / unreadable input")
box(5, 24.2, "Read files in scanner order (parallel chunks)\nGPS time present, finite, varying?")
stop(10.6, 24.2, "exit 3: no / constant /\nnon-finite GPS time")
box(5, 21.8, "Rebuild pulses (return pairing per file)\nusable pulses?")
stop(10.6, 21.8, "exit 4: no usable\npulses")
box(5, 19.4, "Platform check: static scanner (TLS)?\n(lines meet ≤ 5 cm, or ground density ≥ 15×)")
stop(10.6, 19.4, "exit 2: TLS refused\n(v1.1 on branch)")
box(5, 16.8, "Path method per file\n1 pulse lines (multi-return)  2 scan angle (single)\n3 lines + scan-angle gap fill", w=5.4, h=1.9)
stop(10.6, 16.8, "exit 4: single returns\nwithout scan angle,\nno delivered path")
box(5, 13.9, "Quality: median line miss ≤ 0.1 m good · ≤ 0.5 m rough\n> 0.5 m unreliable → retry with scan angle", w=5.6, h=1.6)
ok(1.9, 11.4, "Trajectory CSV\n(3dtrees-trajectory/2)\n+ GeoJSON + report")
box(7.4, 11.4, "Occlusion requested?\npath = delivered (always, if given;\nflag if median miss > 0.1 m) else estimated", w=5.2, h=1.8)
box(7.4, 8.6, "Grid extent + crop (ground: DTM or lowest hits)\nfits 16 GiB budget? else 100 m-buffered tiles", w=5.6, h=1.6)
box(7.4, 6.0, "Route echoes to tile spools (compressed)\ntrace every pulse with OccPy, tiles in parallel", w=5.6, h=1.6)
box(7.4, 3.5, "Classify voxels: 1 filled · 2 empty · 3 occluded · 4 unobserved", w=5.8, h=1.2)
ok(7.4, 1.2, "Occlusion COPC (voxels, classes 1–4)\n+ COG (occluded share per 1 m band, 0.5 m px)")
box(11.65, 8.6, "occlusion error:\nreport section,\nempty products,\npath kept", color="#fde2e1", w=2.4, h=1.9)

for y1, y2 in ((28.2, 27.4), (25.8, 25.0), (23.4, 22.6), (21.0, 20.2), (18.6, 17.8), (15.8, 14.7)):
    arrow(5, y1, 5, y2)
for y in (26.6, 24.2, 21.8, 19.4, 16.8):
    arrow(7.3 if y != 16.8 else 7.7, y, 8.8, y, "no" if y in (24.2, 21.8) else ("yes" if y == 19.4 else None))
arrow(4.2, 13.1, 2.3, 12.1)
arrow(5.8, 13.1, 7.0, 12.3)
arrow(7.4, 10.5, 7.4, 9.4)
arrow(7.4, 7.8, 7.4, 6.8)
arrow(7.4, 5.2, 7.4, 4.1)
arrow(7.4, 2.9, 7.4, 2.0)
arrow(10.25, 8.6, 10.4, 8.6)
ax.text(0.2, 0.2, "Exit codes: 0 done · 2 invalid input / TLS · 3 GPS time unusable · 4 no trajectory possible. "
        "Refusals still write the JSON report.", fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "flow.png", dpi=150)
print("ok", len(table), "stage rows")
