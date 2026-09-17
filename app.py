"""
Sequence Coverage Tool — Streamlit app.

Drop in a BioPharma Finder oligonucleotide-mapping export (CSV or Excel) and get
the sequence-coverage spiral, the linear pile-up map + abundance heat-strip, and
a summary table — computed directly from the numbers. No image, no pixel
tracking, no Firefox. Up to 5 samples can be overlaid.
"""
import io
import os
import streamlit as st
import pandas as pd
import numpy as np

import coverage_core as cc

SAMPLE_CSV = os.path.join(os.path.dirname(__file__), "sample_data", "eGFP_25C.csv")

st.set_page_config(page_title="Sequence Coverage Tool", page_icon="🧬", layout="wide")

st.title("🧬 Sequence Coverage Tool")
st.caption(
    "Coverage maps straight from your BioPharma Finder export — no screenshots, "
    "no pixel tracking, no browser dependency."
)

MAX_SAMPLES = 5

# --------------------------------------------------------------------------
# Cached loading
# --------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_tidy(name: str, data: bytes):
    bio = io.BytesIO(data)
    bio.name = name
    df = cc.read_export(bio)
    tidy = cc.normalise(df)
    return df, tidy


# --------------------------------------------------------------------------
# Sidebar — files + controls
# --------------------------------------------------------------------------
with st.sidebar:
    st.header("1 · Load data")
    files = st.file_uploader(
        "BioPharma Finder export(s) — CSV or Excel",
        type=["csv", "xlsx", "xls", "xlsm"],
        accept_multiple_files=True,
        help="Each file is one sample. Up to 5 samples overlay on the spiral.",
    )
    if files and len(files) > MAX_SAMPLES:
        st.warning(f"Using the first {MAX_SAMPLES} files.")
        files = files[:MAX_SAMPLES]
    has_example = os.path.exists(SAMPLE_CSV)
    use_example = (st.checkbox("Use example data (eGFP 25C)",
                               value=False, disabled=bool(files))
                   if has_example else False)

    st.header("2 · Filters")
    min_conf = st.slider("Min confidence score", 0.0, 100.0, 0.0, 1.0)
    min_area = st.number_input("Min MS area", min_value=0.0, value=0.0, step=1e4, format="%.0f")

# assemble inputs as (name, bytes)
inputs = []
if files:
    inputs = [(f.name, f.getvalue()) for f in files]
elif use_example and os.path.exists(SAMPLE_CSV):
    with open(SAMPLE_CSV, "rb") as fh:
        inputs = [("eGFP_25C.csv", fh.read())]

if not inputs:
    example_hint = " (or tick **Use example data**)" if has_example else ""
    st.info(
        f"👈 Upload a BioPharma Finder oligonucleotide-mapping export to begin"
        f"{example_hint}.\n\n"
        "Expected columns: **Positions** (start-end), **Conf. Score**, **MS Area**, "
        "and an **Oligo** / target-name column. Extra columns are ignored."
    )
    st.stop()

# --------------------------------------------------------------------------
# Per-sample configuration
# --------------------------------------------------------------------------
samples = []
with st.sidebar:
    st.header("3 · Samples")

for fname, fbytes in inputs:
    try:
        raw, tidy = load_tidy(fname, fbytes)
    except Exception as e:  # noqa: BLE001
        st.sidebar.error(f"{fname}: {e}")
        continue

    targets = list(pd.unique(tidy["target_name"]))
    default_name = fname.rsplit(".", 1)[0]
    with st.sidebar.expander(f"⚙️ {default_name}", expanded=len(inputs) == 1):
        sname = st.text_input("Sample name", value=default_name, key=f"n_{fname}")
        target = st.selectbox("Sequence / target", targets, key=f"t_{fname}")
        sub = tidy[tidy["target_name"] == target]
        max_end = int(sub["end"].max()) if len(sub) else 0
        length = st.number_input(
            "Reference length (nt)", min_value=1, value=max(max_end, 1),
            key=f"l_{fname}",
            help="Length of the reference sequence. Defaults to the last covered "
                 "position; set the true length to count the uncovered 3′ end.",
        )
    try:
        samples.append(cc.Sample(sname, tidy, target, length=length,
                                 min_conf=min_conf, min_area=min_area))
    except Exception as e:  # noqa: BLE001
        st.sidebar.error(f"{sname}: {e}")

if not samples:
    st.stop()

# --------------------------------------------------------------------------
# Main — tabs
# --------------------------------------------------------------------------
tab_spiral, tab_linear, tab_summary, tab_data = st.tabs(
    ["🌀 Spiral", "📊 Linear map", "📋 Summary", "🔬 Sequence & data"]
)

with tab_spiral:
    st.subheader("Sequence-coverage spiral")
    if len(samples) > 1:
        st.caption("Samples overlaid — each its own colour, overlaps in red.")
    fig = cc.plot_spiral(samples)
    st.pyplot(fig, use_container_width=True)
    c1, c2, c3 = st.columns(3)
    c1.download_button("⬇️ PNG", cc.fig_bytes(fig, "png", 300),
                       "coverage_spiral.png", "image/png")
    c2.download_button("⬇️ SVG", cc.fig_bytes(fig, "svg"),
                       "coverage_spiral.svg", "image/svg+xml")
    c3.download_button("⬇️ PDF", cc.fig_bytes(fig, "pdf"),
                       "coverage_spiral.pdf", "application/pdf")

with tab_linear:
    st.subheader("Linear coverage map + abundance")
    st.caption("Oligonucleotides stacked by position, coloured by relative "
               "abundance; heat-strip below shows summed abundance per base.")
    for s in samples:
        st.markdown(f"**{s.name}** — {s.target_name}")
        figl = cc.plot_linear(s, title=f"{s.name} — {s.target_name}")
        st.pyplot(figl, use_container_width=True)
        d1, d2, d3 = st.columns(3)
        d1.download_button("⬇️ PNG", cc.fig_bytes(figl, "png", 300),
                           f"{s.name}_linear.png", "image/png", key=f"png_{s.name}")
        d2.download_button("⬇️ SVG", cc.fig_bytes(figl, "svg"),
                           f"{s.name}_linear.svg", "image/svg+xml", key=f"svg_{s.name}")
        d3.download_button("⬇️ PDF", cc.fig_bytes(figl, "pdf"),
                           f"{s.name}_linear.pdf", "application/pdf", key=f"pdf_{s.name}")
        st.divider()

with tab_summary:
    st.subheader("Summary statistics")
    summ = pd.DataFrame([s.summary() for s in samples])
    st.dataframe(summ, use_container_width=True, hide_index=True)
    st.download_button("⬇️ Summary CSV", summ.to_csv(index=False).encode(),
                       "coverage_summary.csv", "text/csv")
    st.caption(
        "MS peaks = detected components (charge states). Unique oligos = distinct "
        "position spans. Coverage % = covered bases ÷ reference length. Mean depth "
        "= average overlapping oligos across covered bases."
    )

with tab_data:
    st.subheader("Reconstructed sequence & per-position data")
    s = samples[0] if len(samples) == 1 else st.selectbox(
        "Sample", samples, format_func=lambda x: x.name)
    seq = s.reconstruct_sequence()
    covered = sum(c != "." for c in seq)
    st.caption(f"{s.target_name}: {len(seq)} nt, {covered} covered "
               f"({100*covered/len(seq):.1f}%). Gaps shown as '.'.")
    st.code("\n".join(seq[i:i+80] for i in range(0, len(seq), 80)) or "(none)")
    st.download_button("⬇️ Sequence (FASTA)",
                       f">{s.name}_{s.target_name}\n{seq}\n".encode(),
                       f"{s.name}_reconstructed.fasta", "text/plain")

    perpos = pd.DataFrame({
        "position": np.arange(1, s.length + 1),
        "base": list(seq),
        "covered": s.cov.astype(int),
        "depth": s.depth.astype(int),
        "rel_abundance": s.con.round(6),
    })
    st.download_button("⬇️ Per-position CSV", perpos.to_csv(index=False).encode(),
                       f"{s.name}_perposition.csv", "text/csv")
    with st.expander("Unique oligonucleotides"):
        st.dataframe(s.species.assign(rel_abundance=s.rel.round(6)),
                     use_container_width=True, hide_index=True)
