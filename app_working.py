"""
app.py  —  qPCR Gene Expression Analyzer
-----------------------------------------
Run locally:   streamlit run app.py
Deploy:        push repo to GitHub, connect at share.streamlit.io
"""

import streamlit as st
import pandas as pd
import io
from process_gene_expression import run_pipeline

# ── Page config ───────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="qPCR Analyzer",
    page_icon="🧬",
    layout="centered",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────

st.markdown("""
<style>
    /* ---- font & background ---- */
    html, body, [class*="css"] {
        font-family: "Inter", "Helvetica Neue", Arial, sans-serif;
    }

    /* ---- hero header ---- */
    .hero {
        background: linear-gradient(135deg, #1a3a5c 0%, #2d6a9f 100%);
        border-radius: 12px;
        padding: 2rem 2.2rem 1.6rem;
        margin-bottom: 1.8rem;
        color: white;
    }
    .hero h1 {
        font-size: 1.9rem;
        font-weight: 700;
        letter-spacing: -0.5px;
        margin: 0 0 0.3rem;
        color: white;
    }
    .hero p {
        font-size: 0.97rem;
        opacity: 0.85;
        margin: 0;
        line-height: 1.55;
    }

    /* ---- section labels ---- */
    .section-label {
        font-size: 0.75rem;
        font-weight: 600;
        letter-spacing: 0.08em;
        text-transform: uppercase;
        color: #6b7280;
        margin: 1.6rem 0 0.5rem;
    }

    /* ---- info card ---- */
    .info-card {
        background: #f0f7ff;
        border-left: 4px solid #2d6a9f;
        border-radius: 0 8px 8px 0;
        padding: 0.85rem 1.1rem;
        font-size: 0.9rem;
        color: #1e3a5f;
        margin-bottom: 1.2rem;
        line-height: 1.55;
    }

    /* ---- results table ---- */
    .results-table th {
        background: #2d6a9f !important;
        color: white !important;
    }

    /* ---- download buttons ---- */
    .stDownloadButton > button {
        border-radius: 8px;
        font-weight: 600;
        padding: 0.5rem 1.2rem;
    }

    /* ---- upload area ---- */
    .uploadedFile {
        border-radius: 8px;
    }

    /* hide Streamlit branding in main area */
    footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ── Hero ──────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="hero">
  <h1>🧬 qPCR Gene Expression Analyzer</h1>
  <p>Upload your raw qPCR spreadsheet (columns A–D).
     Delta and relative expression are calculated automatically,
     and a publication-ready chart is generated for download.</p>
</div>
""", unsafe_allow_html=True)

# ── Format guide ──────────────────────────────────────────────────────────────

with st.expander("📋 Expected spreadsheet format", expanded=False):
    st.markdown("""
**Your `.xlsx` file should have four columns — no header row required:**

| Col A | Col B | Col C | Col D |
|---|---|---|---|
| Organ name *(first row of each group only)* | Sample label | Housekeeping gene (e.g. Cyclo) | Test gene (e.g. CD19-3) |

**Rules:**
- Each organ block is separated by a **blank row**
- The **first named group** in each organ block is treated as the **control**
- Replicate rows have a blank col A and a blank or `cont…` col B
- Columns E and F should be **absent or empty** — they are written by this tool

**Example:**

| A | B | C | D |
|---|---|---|---|
| Liver | NSG-1 | 21.6 | 32.547 |
| | cont | 21.361 | 32.436 |
| | NSG-2 | 20.189 | 31.632 |
| | CD19-LNP-2.5 | 19.761 | 30.937 |
| *(blank row)* | | | |
| Spleen | NSG-1 | 22.361 | 34.913 |
| … | … | … | … |
    """)

# ── File upload ───────────────────────────────────────────────────────────────

st.markdown('<p class="section-label">1 — Upload your data</p>',
            unsafe_allow_html=True)

uploaded = st.file_uploader(
    "Drop your .xlsx file here or click to browse",
    type=["xlsx"],
    label_visibility="collapsed",
)

# ── Optional labels ───────────────────────────────────────────────────────────

st.markdown('<p class="section-label">2 — Confirm sample labels (optional)</p>',
            unsafe_allow_html=True)

col1, col2 = st.columns(2)
with col1:
    hk_name   = st.text_input("Housekeeping gene name", value="Cyclo",
                               help="Used as the column C header label in the output file.")
with col2:
    test_name = st.text_input("Test gene name", value="CD19-3",
                               help="Used as the column D header label in the output file.")

# ── Run ───────────────────────────────────────────────────────────────────────

st.markdown('<p class="section-label">3 — Run analysis</p>',
            unsafe_allow_html=True)

run_btn = st.button("▶  Run Analysis", type="primary", use_container_width=True,
                    disabled=(uploaded is None))

if uploaded is None:
    st.markdown(
        '<div class="info-card">Upload a file above to enable the analysis.</div>',
        unsafe_allow_html=True,
    )

# ── Pipeline ──────────────────────────────────────────────────────────────────

if run_btn and uploaded is not None:
    file_bytes = uploaded.read()

    with st.spinner("Calculating delta and relative expression…"):
        output = run_pipeline(file_bytes)

    if output["errors"]:
        st.error(f"**Error:** {output['errors']}")
        st.stop()

    # ── Results ───────────────────────────────────────────────────────────────

    st.success("Analysis complete!")

    # Summary table — one row per organ × group combination
    st.markdown('<p class="section-label">Results summary</p>',
                unsafe_allow_html=True)

    chart_data = output["chart_data"]
    if chart_data:
        rows = []
        for d in chart_data:
            # Control row
            rows.append({
                "Organ":   d["organ"],
                "Group":   d["ctrl_label"],
                "Type":    "Control",
                "Mean RE": f"{d['ctrl_mean']:.4f}",
            })
            # One row per drug group
            for dg in d["drug_groups"]:
                rows.append({
                    "Organ":   d["organ"],
                    "Group":   dg["label"],
                    "Type":    "Drug",
                    "Mean RE": f"{dg['mean_re']:.4f}",
                })
        df = pd.DataFrame(rows)
        st.dataframe(df, use_container_width=True, hide_index=True)

    # Chart preview
    st.markdown('<p class="section-label">Chart preview</p>',
                unsafe_allow_html=True)

    if output["chart_bytes"]:
        st.image(output["chart_bytes"], use_container_width=True)

    # Downloads
    st.markdown('<p class="section-label">Download results</p>',
                unsafe_allow_html=True)

    dl_col1, dl_col2 = st.columns(2)

    with dl_col1:
        if output["xlsx_bytes"]:
            st.download_button(
                label="⬇  Download Excel file",
                data=output["xlsx_bytes"],
                file_name="gene_expression_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                use_container_width=True,
            )

    with dl_col2:
        if output["chart_bytes"]:
            st.download_button(
                label="⬇  Download chart (PNG)",
                data=output["chart_bytes"],
                file_name="gene_expression_chart.png",
                mime="image/png",
                use_container_width=True,
            )

    # Block details (collapsed)
    with st.expander("🔍 Detected data blocks (debug info)"):
        for blk in output["blocks"]:
            drug_info = "\n".join(
                f"  - `{dg['label']}` — rows {dg['rows']}"
                for dg in blk["drug_groups"]
            )
            st.markdown(
                f"**{blk['organ']}**\n"
                f"- Control: `{blk['ctrl_label']}` — rows {blk['ctrl_rows']}\n"
                f"- Drug groups:\n{drug_info}"
            )

# ── Footer ────────────────────────────────────────────────────────────────────

st.divider()
st.markdown(
    "<p style='text-align:center; color:#9ca3af; font-size:0.82rem;'>"
    "qPCR Analyzer · delta = test gene − HK gene · "
    "relative expression = 2^(avg control delta − delta)"
    "</p>",
    unsafe_allow_html=True,
)
