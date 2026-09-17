"""
coverage_core.py — data-driven sequence-coverage engine for BioPharma Finder
oligonucleotide-mapping exports.

Replaces the old Colab "Sequence Tool", which recovered coverage/abundance by
pixel-tracking a BioPharma Finder (BPF) coverage-map PNG (Firefox-only, brittle).
Here everything is computed directly from the numeric BPF export (CSV/Excel):

  * parse each detected component's Positions (start-end on the reference),
    Conf. Score and MS Area
  * collapse charge states to unique oligonucleotide species
  * per-position coverage + summed relative abundance
  * spiral plot, linear pile-up map + abundance heat-strip, summary stats

No image processing, no browser dependency.
"""
from __future__ import annotations
import re
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap

# --------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------
# BioPharma Finder colours oligo bars by RELATIVE ABUNDANCE in ~10 bands.
# These are a close BPF-style approximation (high abundance -> red, low -> grey);
# replace BPF_BANDS with exact RGB sampled from a real BPF legend for a pixel
# match. Order: HIGHEST abundance first.  Abundance band centres come from the
# original Colab decode table (0.75, 0.35, 0.15, 0.075, ... 0.0005).
BPF_BANDS = [
    ("#8B1A1A", 0.75),    # dark red   - most abundant
    ("#E02323", 0.35),    # red
    ("#F26B1F", 0.15),    # orange
    ("#F7A81B", 0.075),   # amber
    ("#F2E018", 0.035),   # yellow
    ("#8FCf3F", 0.015),   # yellow-green
    ("#28B463", 0.0075),  # green
    ("#17A8A8", 0.0035),  # teal
    ("#2E7DD1", 0.0015),  # blue
    ("#9AA0A6", 0.0005),  # grey       - least abundant
]
# per-sample colours for the spiral overlay (up to 5 samples + combined)
SAMPLE_COLOURS = ["#3955A5", "#F5BA16", "#1EAC4B", "#8E44AD", "#E67E22"]
COMBINED_COLOUR = "#9E1C20"
NOCOV_COLOUR = "#333333"

_BASE_RE = re.compile(r"([ACGUTIacguti])[rd]")
_POS_RE = re.compile(r"\s*(\d+)\s*-\s*(\d+)")
_POS1_RE = re.compile(r"\s*(\d+)\s*$")
_TGT_RE = re.compile(r"\s*(\d+)\s*:")


# --------------------------------------------------------------------------
# Column resolution (robust to BPF export naming variants)
# --------------------------------------------------------------------------
def _find_col(cols, *needles, exclude=None):
    low = {c: c.lower().strip() for c in cols}
    # 1) exact match wins (e.g. "Oligo" beats "Oligo Sequence")
    for need in needles:
        for orig, lc in low.items():
            if lc == need:
                return orig
    # 2) substring match, optionally excluding a token
    for need in needles:
        for orig, lc in low.items():
            if need in lc and (exclude is None or exclude not in lc):
                return orig
    return None


def resolve_columns(df: pd.DataFrame) -> dict:
    cols = list(df.columns)
    return {
        "positions": _find_col(cols, "positions", "position"),
        "conf": _find_col(cols, "conf. score", "conf score", "confidence", "conf"),
        "area": _find_col(cols, "ms area", "area", "abundance", "intensity"),
        "ident": _find_col(cols, "identification", "identity"),
        "oligo_seq": _find_col(cols, "oligo sequence", "sequence"),
        # target name column: prefer "Oligo"/"Protein"/... but never the sequence col
        "target": _find_col(cols, "oligo", "protein", "sequence name", "target",
                            "gene", exclude="sequence"),
    }


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
def _parse_positions(s):
    m = _POS_RE.match(str(s))
    if m:
        return int(m.group(1)), int(m.group(2))
    m = _POS1_RE.match(str(s))
    if m:
        return int(m.group(1)), int(m.group(1))
    return None, None


def _target_prefix(s):
    m = _TGT_RE.match(str(s))
    return m.group(1) if m else None


def decode_bases(seq):
    return "".join(_BASE_RE.findall(str(seq))).upper()


def read_export(file, sheet=0) -> pd.DataFrame:
    """Read a BPF export from a path or file-like (CSV or Excel)."""
    name = getattr(file, "name", str(file)).lower()
    if name.endswith((".xlsx", ".xls", ".xlsm")):
        return pd.read_excel(file, sheet_name=sheet)
    return pd.read_csv(file)


def normalise(df: pd.DataFrame) -> pd.DataFrame:
    """Return a tidy frame: target_name, target_idx, start, end, conf, area, bases."""
    cm_ = resolve_columns(df)
    if not cm_["positions"]:
        raise ValueError("Could not find a 'Positions' column in the export.")
    out = pd.DataFrame()
    pos = df[cm_["positions"]].map(_parse_positions)
    out["start"] = [p[0] for p in pos]
    out["end"] = [p[1] for p in pos]
    out["conf"] = df[cm_["conf"]] if cm_["conf"] else np.nan
    out["area"] = df[cm_["area"]] if cm_["area"] else np.nan
    out["bases"] = df[cm_["oligo_seq"]].map(decode_bases) if cm_["oligo_seq"] else ""
    idx = df[cm_["ident"]].map(_target_prefix) if cm_["ident"] else None
    name = df[cm_["target"]].astype(str) if cm_["target"] else None
    if name is not None:
        out["target_name"] = name
    elif idx is not None:
        out["target_name"] = idx.map(lambda i: f"Seq {i}")
    else:
        out["target_name"] = "sequence"
    out["target_idx"] = idx if idx is not None else "1"
    out = out.dropna(subset=["start", "end"])
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    return out


# --------------------------------------------------------------------------
# Per-sample computation
# --------------------------------------------------------------------------
class Sample:
    def __init__(self, name, tidy, target_name, length=None,
                 min_conf=0.0, min_area=0.0):
        self.name = name
        self.target_name = target_name
        d = tidy[tidy["target_name"] == target_name].copy()
        if min_conf > 0 and d["conf"].notna().any():
            d = d[d["conf"] >= min_conf]
        if min_area > 0 and d["area"].notna().any():
            d = d[d["area"] >= min_area]
        self.n_peaks = len(d)
        # collapse charge states -> unique species
        agg = {"area": ("area", "sum"), "conf": ("conf", "max"),
               "bases": ("bases", "first")}
        self.species = d.groupby(["start", "end"], as_index=False).agg(**agg)
        self.length = int(length) if length else (int(d["end"].max()) if len(d) else 0)
        self._compute()

    def _compute(self):
        L = self.length
        con = np.zeros(L + 1)
        cov = np.zeros(L + 1, dtype=bool)
        depth = np.zeros(L + 1)
        total = self.species["area"].sum() or 1.0
        for _, r in self.species.iterrows():
            s, e = int(r["start"]), min(int(r["end"]), L)
            con[s:e + 1] += r["area"] / total
            cov[s:e + 1] = True
            depth[s:e + 1] += 1
        self.con = con[1:]
        self.cov = cov[1:]
        self.depth = depth[1:]
        self.con_norm = self.con / self.con.max() if self.con.max() > 0 else self.con
        self.rel = self.species["area"] / total
        self.coverage_pct = 100 * self.cov.sum() / L if L else 0.0

    def reconstruct_sequence(self):
        L = self.length
        seq = ["."] * (L + 1)
        for _, r in self.species.iterrows():
            b = r["bases"] or ""
            if len(b) != (int(r["end"]) - int(r["start"]) + 1):
                continue
            for i, ch in enumerate(b):
                p = int(r["start"]) + i
                if 1 <= p <= L and seq[p] == ".":
                    seq[p] = ch
        return "".join(seq[1:])

    def summary(self):
        return {
            "Sample": self.name,
            "Sequence": self.target_name,
            "Length (nt)": self.length,
            "MS peaks": self.n_peaks,
            "Unique oligos": len(self.species),
            "Coverage %": round(self.coverage_pct, 1),
            "Mean depth": round(float(self.depth[self.cov].mean()) if self.cov.any() else 0, 1),
            "Total MS area": f"{self.species['area'].sum():.3e}",
        }


def species_band_colours(sample: Sample):
    """Assign each unique species a BPF band colour by log relative abundance."""
    area = sample.species["area"].to_numpy(dtype=float)
    if len(area) == 0:
        return []
    la = np.log10(area)
    lo, hi = la.min(), la.max()
    span = (hi - lo) or 1.0
    colours = []
    band_cols = [c for c, _ in BPF_BANDS][::-1]  # low -> high
    for v in la:
        idx = int(round((v - lo) / span * (len(band_cols) - 1)))
        colours.append(band_cols[idx])
    return colours


def lane_pack(species: pd.DataFrame):
    lane_end, lanes = [], []
    for _, r in species.sort_values("start").iterrows():
        placed = None
        for ln in range(len(lane_end)):
            if lane_end[ln] < r["start"]:
                placed = ln
                break
        if placed is None:
            placed = len(lane_end)
            lane_end.append(0)
        lane_end[placed] = r["end"]
        lanes.append((r, placed))
    return lanes, len(lane_end)


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def plot_spiral(samples, title=None):
    """Spiral coverage plot (Colab-style). One or more samples overlaid."""
    L = max(s.length for s in samples)
    fig = plt.figure(figsize=(6.6, 8.4))
    t = np.linspace(10 * np.pi, (L * 0.01 + 10) * np.pi, L)

    # combined coverage / abundance across samples
    any_cov = np.zeros(L, dtype=bool)
    con = np.zeros(L)
    per_cov = []
    for s in samples:
        c = np.zeros(L, dtype=bool)
        c[:len(s.cov)] = s.cov
        cc = np.zeros(L)
        cc[:len(s.con_norm)] = s.con_norm
        any_cov |= c
        con = np.maximum(con, cc)
        per_cov.append(c)
    con_norm = con / con.max() if con.max() > 0 else con

    # sample colour per position (single sample -> its colour; multi -> overlap logic)
    if len(samples) == 1:
        colour = np.where(any_cov, SAMPLE_COLOURS[0], NOCOV_COLOUR)
    else:
        colour = np.full(L, NOCOV_COLOUR, dtype=object)
        cov_count = np.sum(per_cov, axis=0)
        for i, c in enumerate(per_cov):
            colour[c] = SAMPLE_COLOURS[i % len(SAMPLE_COLOURS)]
        colour[cov_count > 1] = COMBINED_COLOUR

    con4 = con_norm * 100 + 10
    con5 = (con_norm + 1) / np.max(con_norm + 1)

    # faint radial guide spokes + minor labels
    for i in range(0, 200, 20):
        ss = i
        se = i + ((L // 200) * 200)
        if se >= L:
            se -= 200
        if se <= ss:
            continue
        plt.plot([t[ss] * np.cos(t[ss]), t[se] * np.cos(t[se])],
                 [t[ss] * np.sin(t[ss]), t[se] * np.sin(t[se])],
                 color=[.82, .82, .82], lw=0.9, zorder=1)
        lab = i if i < 100 else i - 100
        if 0 < i < 200:
            plt.text(t[se] * np.cos(t[se]) * 1.06, t[se] * np.sin(t[se]) * 1.06,
                     str(lab), ha="center", va="center", fontsize=8, color="0.45")

    for i in range(L):
        plt.scatter(t[i] * np.cos(t[i]), t[i] * np.sin(t[i]),
                    s=con4[i] * (i / (L / 2) * 0.5 + 0.5),
                    c=colour[i], alpha=con5[i], edgecolors="none", zorder=2)

    plt.scatter(t[0] * np.cos(t[0]), t[0] * np.sin(t[0]), s=60, c="#222", zorder=3)
    plt.text(t[0] * np.cos(t[0]), t[0] * np.sin(t[0]) + 3, "Start",
             fontweight="bold", color="#222", rotation=60)
    div = 100
    for i in range(div, (L // div) * div + 1, div):
        plt.scatter(t[i] * np.cos(t[i]), t[i] * np.sin(t[i]), s=55, c="#222", zorder=3)
        if t[0] * np.cos(t[i]) > 0:
            plt.text(t[i] * np.cos(t[i]), t[i] * np.sin(t[i]) + 3, str(i),
                     fontweight="bold", color="#222", rotation=60)
        else:
            plt.text(t[i] * np.cos(t[i]), t[i] * np.sin(t[i]) - 11, str(i),
                     fontweight="bold", color="#222", ha="right", rotation=60)
    plt.scatter(t[-1] * np.cos(t[-1]), t[-1] * np.sin(t[-1]), s=60, c="#222", zorder=3)
    plt.text(t[-1] * np.cos(t[-1]), t[-1] * np.sin(t[-1]) + 3, f"End ({L})",
             fontweight="bold", color="#222", ha="left", rotation=60)

    plt.axis("off")
    plt.axis("equal")
    plt.tight_layout()
    m = np.sqrt((t[-1] * np.cos(t[-1])) ** 2 + (t[-1] * np.sin(t[-1])) ** 2) * 1.1
    plt.ylim([-1.6 * m, m])
    _rbox = lambda fc: {"facecolor": fc, "edgecolor": "none", "boxstyle": "round"}
    plt.text(0, -m * 1.10, "No coverage", ha="center", va="center", bbox=_rbox("#999"))
    for k, s in enumerate(samples):
        col = SAMPLE_COLOURS[k % len(SAMPLE_COLOURS)]
        plt.text(0, -m * 1.10 - m * 2 / 15 * (k + 1),
                 f"{s.name}: {round(s.coverage_pct)}%", ha="center", va="center",
                 color="white", bbox=_rbox(col))
    if len(samples) > 1:
        plt.text(0, -m * 1.10 - m * 2 / 15 * (len(samples) + 1),
                 f"Combined: {round(100 * any_cov.sum() / L)}%", ha="center",
                 va="center", color="white", bbox=_rbox(COMBINED_COLOUR))
    plt.gca().invert_xaxis()
    if title:
        plt.title(title, fontsize=11)
    return fig


def plot_linear(sample: Sample, title=None):
    """Unwrapped pile-up map + abundance heat-strip (Colab linear output)."""
    from matplotlib import gridspec
    L = sample.length
    lanes, nlanes = lane_pack(sample.species)
    colours = species_band_colours(sample)
    colour_by_idx = {i: colours[i] for i in range(len(colours))}
    # map species row -> colour via positional index
    sp = sample.species.reset_index(drop=True)
    idx_of = {(int(r.start), int(r.end)): i for i, r in sp.iterrows()}

    fig = plt.figure(figsize=(18, max(2.6, 0.13 * nlanes + 1.4)))
    gs = gridspec.GridSpec(2, 1, height_ratios=[max(4, nlanes * 0.5), 1], hspace=0.18)
    ax = fig.add_subplot(gs[0])
    ax.set_facecolor((0.922, 0.933, 0.964))
    for (r, ln) in lanes:
        col = colour_by_idx.get(idx_of[(int(r["start"]), int(r["end"]))], "#888")
        ax.add_patch(Rectangle((r["start"], -(ln + 1)), r["end"] - r["start"] + 1, 0.8,
                               facecolor=col, edgecolor="white", linewidth=0.2))
    ax.set_xlim(0, L + 1)
    ax.set_ylim(-(nlanes + 1), 1)
    ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False,
                   labelsize=8, length=2, colors="#555")
    ax.set_xticks(range(100, L + 1, 100))
    ax.set_title(title or f"{sample.name} — oligonucleotide coverage",
                 loc="left", fontsize=11, pad=18)

    axh = fig.add_subplot(gs[1])
    blues = LinearSegmentedColormap.from_list("b", ["#ffffff", SAMPLE_COLOURS[0]])
    axh.imshow(sample.con_norm.reshape(1, -1), aspect="auto", cmap=blues,
               extent=[1, L, 0, 1], vmin=0, vmax=1)
    axh.set_yticks([])
    axh.set_xlim(0, L + 1)
    axh.set_xticks(range(100, L + 1, 100))
    axh.tick_params(labelsize=8, length=2, colors="#555")
    for sp_ in axh.spines.values():
        sp_.set_visible(False)
    axh.set_ylabel("abundance", fontsize=8, rotation=0, ha="right", va="center")
    return fig


def fig_bytes(fig, fmt="png", dpi=200):
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()
