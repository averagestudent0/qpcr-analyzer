import io
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

COL_ORGAN = 1   # A
COL_LABEL = 2   # B
COL_HK    = 3   # C
COL_TEST  = 4   # D
COL_DELTA = 5   # E
COL_REL   = 6   # F

CTRL_HEX = "D9D9D9"
DRUG_HEX = "F4A99A"
HDR_HEX  = "2F4F8F"
CTRL_MPL = "#C8C8C8"
DRUG_MPL = "#E04030"


def cl(n):
    return get_column_letter(n)

def is_numeric(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False

def is_continuation(label):
    if label is None:
        return True
    s = str(label).strip().lower()
    return s == "" or s.startswith("cont")


# ── Parser ────────────────────────────────────────────────────────────────────

def parse_blocks(ws) -> list:
    """
    Return one dict per organ block:
    {
        organ, hk_name, test_name,
        header_row (the text/label row just above first data row, or None),
        first_data_row, last_data_row,
        ctrl_rows, ctrl_label,
        drug_rows, drug_labels,
    }
    """
    blocks     = []
    max_row    = ws.max_row

    # Identify every row that has numeric values in both C and D
    numeric_rows = set()
    for r in range(1, max_row + 1):
        if is_numeric(ws.cell(row=r, column=COL_HK).value) and \
           is_numeric(ws.cell(row=r, column=COL_TEST).value):
            numeric_rows.add(r)

    visited = set()
    r = 1
    while r <= max_row:
        if r not in numeric_rows or r in visited:
            r += 1
            continue

        # ── find the contiguous block of numeric rows starting at r ──
        block_start = r
        block_end   = r
        while block_end + 1 in numeric_rows:
            block_end += 1

        for rr in range(block_start, block_end + 1):
            visited.add(rr)

        # ── look for organ name / column headers in the row(s) above ──
        organ     = None
        hk_name   = None
        test_name = None

        for look in range(block_start - 1, max(0, block_start - 4), -1):
            a = ws.cell(row=look, column=COL_ORGAN).value
            c = ws.cell(row=look, column=COL_HK).value
            d = ws.cell(row=look, column=COL_TEST).value
            if a is not None and str(a).strip():
                organ = str(a).strip()
            if isinstance(c, str) and c.strip():
                hk_name = c.strip()
            if isinstance(d, str) and d.strip():
                test_name = d.strip()
            if organ:
                break

        # Also check col A within the data rows themselves
        if organ is None:
            for rr in range(block_start, block_end + 1):
                a = ws.cell(row=rr, column=COL_ORGAN).value
                if a is not None and str(a).strip():
                    organ = str(a).strip()
                    break

        if organ is None:
            organ = f"Block_{block_start}"

        # ── group data rows by sample ──
        groups      = []
        current_grp = None

        for rr in range(block_start, block_end + 1):
            bv = ws.cell(row=rr, column=COL_LABEL).value
            if is_continuation(bv):
                if current_grp is not None:
                    current_grp["rows"].append(rr)
            else:
                current_grp = {"name": str(bv).strip(), "rows": [rr]}
                groups.append(current_grp)

        if not groups:
            r = block_end + 1
            continue

        ctrl_rows   = groups[0]["rows"]
        ctrl_label  = groups[0]["name"]
        drug_rows   = []
        drug_labels = []
        for g in groups[1:]:
            drug_rows.extend(g["rows"])
            drug_labels.append(g["name"])

        blocks.append({
            "organ":          organ,
            "hk_name":        hk_name   or "HK Gene",
            "test_name":      test_name or "Test Gene",
            "first_data_row": block_start,
            "last_data_row":  block_end,
            "ctrl_rows":      ctrl_rows,
            "ctrl_label":     ctrl_label,
            "drug_rows":      drug_rows,
            "drug_labels":    drug_labels,
        })

        r = block_end + 1

    return blocks


# ── Numerical computation ─────────────────────────────────────────────────────

def compute_values(ws, blocks) -> list:
    """
    Compute delta and relative expression using raw float values.
    Augments each block dict in-place and returns the list.
    """
    for blk in blocks:
        all_rows = blk["ctrl_rows"] + blk["drug_rows"]

        deltas = {}
        for row in all_rows:
            hk   = float(ws.cell(row=row, column=COL_HK).value)
            test = float(ws.cell(row=row, column=COL_TEST).value)
            deltas[row] = test - hk

        ctrl_deltas   = [deltas[r] for r in blk["ctrl_rows"]]
        avg_ctrl      = sum(ctrl_deltas) / len(ctrl_deltas)
        rel           = {r: 2 ** (avg_ctrl - deltas[r]) for r in all_rows}

        ctrl_rel = [rel[r] for r in blk["ctrl_rows"]]
        drug_rel = [rel[r] for r in blk["drug_rows"]]

        blk["deltas"]         = deltas
        blk["avg_ctrl_delta"] = avg_ctrl
        blk["ctrl_rel"]       = ctrl_rel
        blk["drug_rel"]       = drug_rel
        blk["ctrl_mean_re"]   = float(np.mean(ctrl_rel)) if ctrl_rel else 0.0
        blk["drug_mean_re"]   = float(np.mean(drug_rel)) if drug_rel else 0.0

    return blocks


# ── Output workbook builder ───────────────────────────────────────────────────

def build_output_workbook(ws_source, blocks) -> bytes:
    """
    Clone the source worksheet and insert header rows + formulas for cols E & F.
    Returns raw xlsx bytes of the completed workbook.
    """
    # Save source to bytes and reload fresh
    src_buf = io.BytesIO()
    ws_source.parent.save(src_buf)
    wb_out = load_workbook(io.BytesIO(src_buf.getvalue()))
    ws_out = wb_out.active

    # Insert header rows from bottom to top so earlier row indices stay valid
    sorted_blocks = sorted(blocks, key=lambda b: b["first_data_row"], reverse=True)
    offset_map    = {}   # organ → cumulative offset at time of processing

    # We insert top-down so process top-to-bottom with a running offset
    sorted_asc    = sorted(blocks, key=lambda b: b["first_data_row"])
    cumulative    = 0

    for blk in sorted_asc:
        insert_at = blk["first_data_row"] + cumulative
        ws_out.insert_rows(insert_at)

        # Write column labels into header row
        ws_out.cell(row=insert_at, column=COL_HK).value    = blk["hk_name"]
        ws_out.cell(row=insert_at, column=COL_TEST).value  = blk["test_name"]
        ws_out.cell(row=insert_at, column=COL_DELTA).value = "delta"
        ws_out.cell(row=insert_at, column=COL_REL).value   = "Avg Ctrl Δ / Rel Expr"

        blk["_hdr"]       = insert_at
        blk["_ctrl_out"]  = [r + cumulative + 1 for r in blk["ctrl_rows"]]
        blk["_drug_out"]  = [r + cumulative + 1 for r in blk["drug_rows"]]
        cumulative       += 1

    # Write formulas
    for blk in sorted_asc:
        h       = blk["_hdr"]
        all_out = blk["_ctrl_out"] + blk["_drug_out"]

        for row in all_out:
            ws_out.cell(row=row, column=COL_DELTA).value = (
                f"={cl(COL_TEST)}{row}-{cl(COL_HK)}{row}"
            )

        ctrl_refs = ",".join(f"{cl(COL_DELTA)}{r}" for r in blk["_ctrl_out"])
        ws_out.cell(row=h, column=COL_REL).value = f"=AVERAGE({ctrl_refs})"

        for row in all_out:
            ws_out.cell(row=row, column=COL_REL).value = (
                f"=2^({cl(COL_REL)}{h}-{cl(COL_DELTA)}{row})"
            )

    # Add legend row at very top
    ws_out.insert_rows(1)
    ws_out["A1"] = "Legend:"
    ws_out["B1"] = "Gray = Control   |   Red = Drug-treated"
    ws_out["A1"].font = Font(bold=True, name="Arial", size=10)
    ws_out["B1"].font = Font(italic=True, name="Arial", size=10)

    # Style (all row refs shift +1 for legend row)
    _style_output(ws_out, sorted_asc, legend_offset=1)

    out_buf = io.BytesIO()
    wb_out.save(out_buf)
    return out_buf.getvalue()


# ── Styling ───────────────────────────────────────────────────────────────────

def _fill(hex_color):
    return PatternFill("solid", start_color=hex_color, end_color=hex_color)

def _border():
    t = Side(style="thin", color="CCCCCC")
    return Border(left=t, right=t, top=t, bottom=t)

_CENTER = Alignment(horizontal="center", vertical="center")

def _style_output(ws, blocks, legend_offset=1):
    ws.column_dimensions[cl(COL_ORGAN)].width = 10
    ws.column_dimensions[cl(COL_LABEL)].width = 22
    ws.column_dimensions[cl(COL_HK)].width    = 14
    ws.column_dimensions[cl(COL_TEST)].width   = 14
    ws.column_dimensions[cl(COL_DELTA)].width  = 12
    ws.column_dimensions[cl(COL_REL)].width    = 22

    hdr_fill  = _fill(HDR_HEX)
    ctrl_fill = _fill(CTRL_HEX)
    drug_fill = _fill(DRUG_HEX)
    border    = _border()
    hdr_font  = Font(bold=True, color="FFFFFF", name="Arial", size=10)
    std_font  = Font(name="Arial", size=10)

    for blk in blocks:
        h         = blk["_hdr"]  + legend_offset
        ctrl_rows = [r + legend_offset for r in blk["_ctrl_out"]]
        drug_rows = [r + legend_offset for r in blk["_drug_out"]]

        for col in range(1, COL_REL + 1):
            cell = ws.cell(row=h, column=col)
            cell.font = hdr_font; cell.fill = hdr_fill
            cell.alignment = _CENTER; cell.border = border

        for row in ctrl_rows:
            for col in range(1, COL_REL + 1):
                cell = ws.cell(row=row, column=col)
                cell.font = std_font; cell.fill = ctrl_fill
                cell.alignment = _CENTER; cell.border = border

        for row in drug_rows:
            for col in range(1, COL_REL + 1):
                cell = ws.cell(row=row, column=col)
                cell.font = std_font; cell.fill = drug_fill
                cell.alignment = _CENTER; cell.border = border


# ── Chart renderer ────────────────────────────────────────────────────────────

def render_chart(blocks: list) -> bytes:
    organs     = [b["organ"]        for b in blocks]
    ctrl_means = [b["ctrl_mean_re"] for b in blocks]
    drug_means = [b["drug_mean_re"] for b in blocks]
    n          = len(organs)

    ctrl_label = blocks[0]["ctrl_label"]             if blocks else "Control"
    drug_label = ", ".join(blocks[0]["drug_labels"]) if blocks else "Treatment"

    group_spacing = 1.4
    bar_width     = 0.22
    half_span     = bar_width / 2 + 0.05
    x             = np.arange(n) * group_spacing

    fig, ax = plt.subplots(figsize=(max(5, 2.4 * n + 1.8), 5.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    ax.bar(x - half_span, ctrl_means, width=bar_width, color=CTRL_MPL, zorder=3)
    ax.bar(x + half_span, drug_means, width=bar_width, color=DRUG_MPL, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(organs, fontsize=13)
    ax.set_ylabel("Relative Expression", fontsize=12)
    ax.set_ylim(bottom=0)
    ax.yaxis.grid(True, color="#E8E8E8", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#AAAAAA")
    ax.spines["bottom"].set_color("#AAAAAA")
    ax.tick_params(axis="both", which="both", length=0, labelsize=11)
    ax.set_xlim(-group_spacing * 0.55,
                (n - 1) * group_spacing + group_spacing * 0.55)

    ctrl_patch = mpatches.Patch(color=CTRL_MPL, label=ctrl_label)
    drug_patch = mpatches.Patch(color=DRUG_MPL, label=drug_label)
    ax.legend(handles=[ctrl_patch, drug_patch], loc="upper center",
              bbox_to_anchor=(0.5, 1.13), ncol=2, frameon=False,
              fontsize=11, handlelength=1.0, handleheight=1.0)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=180, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close()
    buf.seek(0)
    return buf.read()


# ── Top-level pipeline ────────────────────────────────────────────────────────

def run_pipeline(file_bytes: bytes) -> dict:
    result = {
        "xlsx_bytes": None, "chart_bytes": None,
        "chart_data": [],   "blocks": [],   "errors": "",
    }

    try:
        wb = load_workbook(io.BytesIO(file_bytes))
        ws = wb.active
    except Exception as e:
        result["errors"] = f"Could not open file: {e}"
        return result

    blocks = parse_blocks(ws)
    if not blocks:
        result["errors"] = (
            "No data blocks detected. Make sure columns C and D contain "
            "numeric values and each organ group is separated by a blank row."
        )
        return result

    try:
        blocks = compute_values(ws, blocks)
    except Exception as e:
        result["errors"] = f"Calculation error: {e}"
        return result

    result["blocks"] = blocks
    result["chart_data"] = [
        {
            "organ":       b["organ"],
            "ctrl_mean":   b["ctrl_mean_re"],
            "drug_mean":   b["drug_mean_re"],
            "ctrl_label":  b["ctrl_label"],
            "drug_labels": b["drug_labels"],
        }
        for b in blocks
    ]

    try:
        result["xlsx_bytes"] = build_output_workbook(ws, blocks)
    except Exception as e:
        result["errors"] = f"Excel output error: {e}"
        return result

    try:
        result["chart_bytes"] = render_chart(blocks)
    except Exception as e:
        result["errors"] = f"Chart rendering failed: {e}"
        return result

    return result
