# 🧬 Sequence Coverage Tool

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://sequence-coverage-tool.streamlit.app/)

**🔗 Live app: [sequence-coverage-tool.streamlit.app](https://sequence-coverage-tool.streamlit.app/)**

A [Streamlit](https://streamlit.io/) app that turns a **BioPharma Finder**
oligonucleotide-mapping export into sequence-coverage visualisations —
**computed directly from the numbers**, with no image processing, no pixel
tracking, and no Firefox dependency.

This replaces the old Google Colab "Sequence Tool", which recovered coverage and
abundance by reading a BioPharma Finder coverage-map **PNG** pixel-by-pixel (only
rendered correctly in Firefox, and brittle). Everything here is derived from the
BioPharma Finder CSV/Excel export instead, so results are exact and reproducible.

## What it does

Upload one or more BioPharma Finder exports (CSV or Excel) and get:

- **🌀 Spiral map** — the sequence wound as a spiral, each base coloured by
  coverage and sized/faded by abundance. Up to **5 samples** overlay (each its
  own colour, overlaps in red).
- **📊 Linear map** — oligonucleotides stacked by position and coloured by
  relative abundance (BioPharma Finder-style bands), with an abundance
  heat-strip beneath.
- **📋 Summary** — MS peaks, unique oligos, coverage %, mean depth, total MS
  area per sample.
- **🔬 Sequence & data** — the reference sequence reconstructed from the mapped
  oligos, plus per-position and per-oligo tables. Everything is downloadable
  (PNG / SVG / PDF / CSV / FASTA).

## Input

Any BioPharma Finder component export with these columns (extra columns ignored,
names matched flexibly):

| Needs | Column (typical) |
|-------|------------------|
| Position span | `Positions` (e.g. `372-376`) |
| Confidence | `Conf. Score` |
| Abundance | `MS Area` |
| Target name | `Oligo` (e.g. `eGFP`, `eGFP_Randomized_1`) |
| Oligo sequence *(optional)* | `Oligo Sequence` (used to reconstruct the reference) |

Multiple charge states of the same oligo are collapsed to one species by
position span.

> **Example data:** drop a BioPharma Finder export at `sample_data/eGFP_25C.csv`
> to enable the built-in **Use example data** toggle. (Not bundled in this repo.)

## Run it

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy / update

This repo is deployed on [Streamlit Community Cloud](https://share.streamlit.io/)
from the `main` branch (`app.py` as the entry point), live at
**[sequence-coverage-tool.streamlit.app](https://sequence-coverage-tool.streamlit.app/)**.
Any push to `main` **auto-redeploys** the live app within a minute or two — no
manual redeploy step.

## Notes

- **Reference length**: coverage % needs the true sequence length. The app
  defaults to the last covered position; set the real length (e.g. 882 for the
  eGFP construct) to count the uncovered 3′ end.
- **Colours**: bars use a BioPharma Finder-style abundance palette in
  `coverage_core.py` (`BPF_BANDS`). Swap in exact RGB sampled from a BioPharma
  Finder legend for a pixel-perfect match.
