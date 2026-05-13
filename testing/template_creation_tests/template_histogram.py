"""
compare_template_pccs.py

For every subject that has Daniel, Mira, and Deep templates, compute the
per-bin pairwise PCC between each pair of sources and plot the distribution
of (1 - PCC) -- i.e. distance from perfect agreement, so exact matches sit
at 0 like the SSD plots in compare_wave_bounds.py.

Emits, for each of the three pairs (Daniel-Mira, Daniel-Deep, Mira-Deep):
  * <pair>_per_subject_histograms.svg : 5-column per-subject grid, log-y,
    bucket 0 reserved for exact-zero values, plus a zero-diff reference panel
  * <pair>_pooled_histogram.svg       : one chart pooling all bins across
    every subject

Template reading and the PCC routine are imported from visualize_templates.py
so this script stays a thin layer of "iterate + score + plot."
"""

import math as _math
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Reuse readers, paths, and the alignment-aware correlation from the
# existing visualizer. corr_pair_aligned is what plot_chunk already uses
# to compute the three numbers printed in each panel title, so the scores
# in these histograms match those titles exactly.
# ---------------------------------------------------------------------------
from visualize_templates import (
    DANIEL_SR,
    DEEP_SR_DEFAULT,
    MIRA_SR,
    corr_pair_aligned,
    daniel_path,
    deep_path,
    find_matching_files,
    mira_path,
    normalize,
    read_daniel_template_mat,
    read_deep_template_mat,
    read_mira_template_bin,
)

OUTPUT_DIR = Path(r"D:\USERS\MiraWelner\QTVI\testing\template_creation_tests")


# ============================================================================
# Per-subject scoring
# ============================================================================


def _vec(t, key):
    if t is None:
        return np.array([], dtype=np.float64)
    return np.asarray(t.get(key, []), dtype=np.float64)


def score_subject(daniel_list, mira_list, deep_list, deep_sr):
    """Return three lists (Daniel-Mira, Daniel-Deep, Mira-Deep) of per-bin
    1-PCC values. Bins where a PCC is NaN (insufficient data or zero-variance
    template, or bin flagged as bad on the Mira side) are skipped.
    """
    dm, dp, mp = [], [], []
    n = max(len(daniel_list), len(mira_list), len(deep_list))

    for i in range(n):
        d = daniel_list[i] if i < len(daniel_list) else None
        m = mira_list[i] if i < len(mira_list) else None
        p = deep_list[i] if i < len(deep_list) else None

        # Skip the same way the visualizer does: Mira-flagged bad bins
        # produce meaningless templates and would just pollute the
        # distribution.
        if m is not None and m.get("bad_segment", False):
            continue

        d_vec = _vec(d, "ecgTemplate")
        m_vec = _vec(m, "ch1_raw_template")
        p_vec = _vec(p, "avgcyc")

        d_n = normalize(d_vec)
        m_n = normalize(m_vec)
        p_n = normalize(p_vec)

        c_dm = corr_pair_aligned(d_n, m_n, sr_a=DANIEL_SR, sr_b=MIRA_SR)
        c_dp = corr_pair_aligned(d_n, p_n, sr_a=DANIEL_SR, sr_b=deep_sr)
        c_mp = corr_pair_aligned(m_n, p_n, sr_a=MIRA_SR, sr_b=deep_sr)

        # 1 - PCC, clipped at 0 so floating-point fuzz on perfect matches
        # doesn't push a 1.0000000002 correlation slightly below zero on
        # the distance axis.
        if not np.isnan(c_dm):
            dm.append(max(0.0, 1.0 - float(c_dm)))
        if not np.isnan(c_dp):
            dp.append(max(0.0, 1.0 - float(c_dp)))
        if not np.isnan(c_mp):
            mp.append(max(0.0, 1.0 - float(c_mp)))

    return dm, dp, mp


# ============================================================================
# SVG helpers (matches the style of compare_wave_bounds.py)
# ============================================================================


def _make_log_y_axis(lines, cx, cy, chart_h, max_count):
    log_max = _math.log10(max_count + 1) * 1.3 if max_count > 0 else 1.0

    def bar_h(count):
        if count <= 0:
            return 0.0
        return chart_h * _math.log10(count + 1) / log_max

    tick_vals = [0]
    v = 1
    while v <= max_count + 1:
        tick_vals.append(v)
        v *= 10
    if max_count not in tick_vals:
        tick_vals.append(max_count)

    for tv in tick_vals:
        yp = cy + chart_h - bar_h(tv)
        if yp < cy - 2:
            continue
        is_max = tv == max_count and max_count not in [0, 1, 10, 100, 1000, 10000]
        lines.append(
            f'<line x1="{cx - 3}" y1="{yp:.1f}" x2="{cx + 1}" y2="{yp:.1f}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        bold_attrs = ' fill="#e74c3c" font-weight="bold"' if is_max else ""
        lines.append(
            f'<text x="{cx - 5}" y="{yp + 4:.1f}" '
            f'text-anchor="end" class="tick"{bold_attrs}>{tv}</text>'
        )
        if tv > 0:
            lines.append(
                f'<line x1="{cx + 1}" y1="{yp:.1f}" x2="{cx + 9999}" y2="{yp:.1f}" '
                f'stroke="#eee" stroke-width="0.5" stroke-dasharray="2,2"/>'
            )
    return bar_h


def _bucket_counts(vals, max_val, num_buckets):
    """Bucket 0 = exact zero. Remaining buckets uniformly tile (0, max_val]."""
    counts = np.zeros(num_buckets, dtype=int)
    if max_val <= 0:
        # All zeros (or empty): everything in bucket 0.
        counts[0] = sum(1 for v in vals if v == 0.0)
        return counts
    n_nonzero_bins = num_buckets - 1
    bucket_width = max_val / n_nonzero_bins
    for v in vals:
        if v == 0.0:
            counts[0] += 1
        else:
            bi = 1 + int((v - 1e-300) / bucket_width)
            if bi >= num_buckets:
                bi = num_buckets - 1
            counts[bi] += 1
    return counts


def _append_zero_ref_panel(
    lines,
    ox,
    oy,
    cell_w,
    cell_h,
    pad_l,
    pad_r,
    pad_t,
    pad_b,
    n_total,
    ref_max_val,
    x_label,
    num_buckets=50,
    fill_color="#2ecc71",
):
    chart_w = cell_w - pad_l - pad_r
    chart_h = cell_h - pad_t - pad_b
    cx = ox + pad_l
    cy = oy + pad_t

    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
        f'text-anchor="middle" class="subtitle" fill="#27ae60">'
        f"ZERO-DIFF REFERENCE (n={n_total})</text>"
    )
    lines.append(
        f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
        f'text-anchor="middle" class="stats">'
        f"all 1-PCC=0 | x-range=[0, {ref_max_val:.6f}]</text>"
    )

    lines.append(
        f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<line x1="{cx}" y1="{cy + chart_h}" '
        f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<text x="{ox + 18}" y="{cy + chart_h // 2}" '
        f'text-anchor="middle" class="y-axis-label" '
        f'transform="rotate(-90 {ox + 10},{cy + chart_h // 2})">'
        f"# bins</text>"
    )

    y_axis_max = n_total if n_total > 0 else 1
    for t in range(4):
        y_val = int(round(y_axis_max * t / 3.0))
        y_pos = cy + chart_h - (y_val / y_axis_max) * chart_h
        lines.append(
            f'<text x="{cx - 3}" y="{y_pos + 3}" '
            f'text-anchor="end" class="tick">{y_val}</text>'
        )
        if t > 0:
            lines.append(
                f'<line x1="{cx + 1}" y1="{y_pos}" '
                f'x2="{cx + chart_w}" y2="{y_pos}" '
                f'stroke="#eee" stroke-width="0.5"/>'
            )

    bar_w = chart_w / num_buckets
    bx = cx
    bar_h = chart_h
    by = cy
    lines.append(
        f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
        f'width="{bar_w - 1:.1f}" height="{bar_h:.1f}" '
        f'fill="{fill_color}" fill-opacity="0.8">'
        f"<title>[0, 0]: {n_total}</title></rect>"
    )
    lines.append(
        f'<text x="{bx + bar_w / 2:.1f}" y="{cy - 2}" '
        f'text-anchor="middle" class="clipped" fill="{fill_color}">{n_total}</text>'
    )

    for t in range(3):
        val = ref_max_val * t / 2.0
        x_pos = cx + chart_w * t / 2.0
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="x-tick">{val:.6f}</text>'
        )

    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
        f'text-anchor="middle" class="axis-label">{x_label}</text>'
    )


# ============================================================================
# Per-subject grid (5 cols + zero-ref panel)
# ============================================================================


def write_per_subject_pcc_grid(
    path,
    per_subject_vals,
    chart_title,
    bar_color="#1a5fa8",
    ref_color="#2ecc71",
    num_buckets=50,
    x_label="1 - PCC",
):
    """`per_subject_vals` is a list of (subject_id, list_of_1mPCC) tuples,
    one entry per subject. One panel per subject, plus one zero-diff
    reference panel in the next free slot."""
    if not per_subject_vals:
        return

    num_files = len(per_subject_vals)
    n_cols = 5
    n_rows = (num_files + 1 + n_cols - 1) // n_cols

    cell_w, cell_h = 260, 275
    pad_l, pad_r, pad_t, pad_b = 50, 25, 48, 42
    chart_w = cell_w - pad_l - pad_r
    chart_h = cell_h - pad_t - pad_b

    svg_w = n_cols * cell_w
    svg_h = n_rows * cell_h + 55

    lines = []
    lines.append(
        f'<svg width="{svg_w}" height="{svg_h}" xmlns="http://www.w3.org/2000/svg">'
    )
    lines.append(
        "<style>\n"
        "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        "  .title { font-size: 16px; font-weight: bold; }\n"
        "  .subtitle { font-size: 12px; font-weight: bold; }\n"
        "  .stats { font-size: 12px; fill: #555; }\n"
        "  .tick { font-size: 12px; }\n"
        "  .axis-label { font-size: 12px; }\n"
        "  .y-axis-label { font-size: 12px; }\n"
        "  .x-tick { font-size: 12px; }\n"
        "  .clipped { font-size: 12px; fill: #e74c3c; font-weight: bold; }\n"
        "  .empty { font-size: 12px; fill: #999; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')
    lines.append(
        f'<text x="{svg_w // 2}" y="30" text-anchor="middle" class="title">'
        f"{chart_title} - all {num_files} subjects</text>"
    )

    all_max_vals = []
    all_n = []

    for f_idx, (subject_id, raw_vals) in enumerate(per_subject_vals):
        col = f_idx % n_cols
        row = f_idx // n_cols
        ox = col * cell_w
        oy = row * cell_h + 50

        vals = [round(v, 6) for v in raw_vals]
        n_total = len(vals)
        zero_count = sum(1 for v in vals if v == 0.0)
        zero_pct = (zero_count / n_total * 100.0) if n_total else 0.0

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 14}" '
            f'text-anchor="middle" class="subtitle">'
            f"{subject_id} ({zero_pct:.1f}% 0.000000)</text>"
        )

        if not vals:
            lines.append(
                f'<text x="{ox + cell_w // 2}" y="{oy + cell_h // 2}" '
                f'text-anchor="middle" class="empty">(no bins)</text>'
            )
            continue

        arr = np.array(vals)
        mn = float(np.mean(arr))
        md = float(np.median(arr))

        nonzero = arr[arr > 0]
        max_val = float(nonzero.max()) if len(nonzero) > 0 else 0.0

        all_max_vals.append(max_val if max_val > 0 else 1.0)
        all_n.append(n_total)

        lines.append(
            f'<text x="{ox + cell_w // 2}" y="{oy + 26}" '
            f'text-anchor="middle" class="stats">'
            f"n={n_total}  mean={mn:.6f}  med={md:.6f}</text>"
        )

        counts = _bucket_counts(vals, max_val, num_buckets)
        max_count = int(counts.max()) if counts.max() > 0 else 1

        cx = ox + pad_l
        cy = oy + pad_t
        bar_w = chart_w / num_buckets

        lines.append(
            f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<line x1="{cx}" y1="{cy + chart_h}" '
            f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
            f'stroke="#333" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{ox + 18}" y="{cy + chart_h // 2}" '
            f'text-anchor="middle" class="y-axis-label" '
            f'transform="rotate(-90 {ox + 10},{cy + chart_h // 2})">'
            f"# bins (log)</text>"
        )

        bar_h_fn = _make_log_y_axis(lines, cx, cy, chart_h, max_count)

        # If max_val is 0 (every bin matched perfectly), only bucket 0 will
        # have anything in it; bucket_width below is irrelevant since the
        # loop bails on c == 0 anyway.
        n_nonzero_bins = num_buckets - 1
        bucket_width = (
            (max_val / n_nonzero_bins) if (max_val > 0 and n_nonzero_bins > 0) else 0.0
        )

        for b_idx in range(num_buckets):
            c = counts[b_idx]
            if c == 0:
                continue
            bh = bar_h_fn(c)
            bx = cx + b_idx * bar_w
            by = cy + chart_h - bh
            if b_idx == 0:
                tip = f"[0, 0]: {c}"
            else:
                rs = (b_idx - 1) * bucket_width
                re = b_idx * bucket_width
                tip = f"({rs:.6f}, {re:.6f}]: {c}"
            lines.append(
                f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
                f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
                f'fill="{bar_color}" fill-opacity="0.92">'
                f"<title>{tip}</title></rect>"
            )

        # X-axis ticks: 0, max/2, max -- only meaningful when max_val > 0.
        for t in range(3):
            val = (max_val if max_val > 0 else 1.0) * t / 2.0
            x_pos = cx + chart_w * t / 2.0
            lines.append(
                f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
                f'text-anchor="middle" class="x-tick">{val:.6f}</text>'
            )

        lines.append(
            f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 30}" '
            f'text-anchor="middle" class="axis-label">{x_label}</text>'
        )

    if all_max_vals:
        global_max = float(max(all_max_vals))
        mean_n = int(round(sum(all_n) / len(all_n))) if all_n else 0
        ref_idx = num_files
        ref_col = ref_idx % n_cols
        ref_row = ref_idx // n_cols
        ox_ref = ref_col * cell_w
        oy_ref = ref_row * cell_h + 50
        _append_zero_ref_panel(
            lines,
            ox=ox_ref,
            oy=oy_ref,
            cell_w=cell_w,
            cell_h=cell_h,
            pad_l=pad_l,
            pad_r=pad_r,
            pad_t=pad_t,
            pad_b=pad_b,
            n_total=mean_n,
            ref_max_val=global_max,
            x_label=x_label,
            num_buckets=num_buckets,
            fill_color=ref_color,
        )

    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Pooled histogram (one big chart, all subjects combined)
# ============================================================================


def write_pooled_pcc_histogram(
    path,
    per_subject_vals,
    chart_title,
    num_buckets=80,
    bar_color="#3498db",
    x_label="1 - PCC",
):
    all_vals = [round(v, 6) for _, vals in per_subject_vals for v in vals]
    if not all_vals:
        return

    arr = np.array(all_vals)
    mn = float(np.mean(arr))
    md = float(np.median(arr))
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    zero_count = sum(1 for v in all_vals if v == 0.0)
    zero_pct = zero_count / len(all_vals) * 100.0

    svg_w, svg_h = 900, 500
    pad_l, pad_r, pad_t, pad_b = 70, 40, 70, 60
    chart_w = svg_w - pad_l - pad_r
    chart_h = svg_h - pad_t - pad_b

    nonzero = arr[arr > 0]
    max_val = float(nonzero.max()) if len(nonzero) > 0 else 0.0

    counts = _bucket_counts(all_vals, max_val, num_buckets)
    max_count = int(counts.max()) if counts.max() > 0 else 1
    bar_w = chart_w / num_buckets
    cx, cy = pad_l, pad_t

    n_subjects = len(per_subject_vals)

    lines = []
    lines.append(
        f'<svg width="{svg_w}" height="{svg_h}" xmlns="http://www.w3.org/2000/svg">'
    )
    lines.append(
        "<style>\n"
        "  text { font-family: Consolas, 'Courier New', monospace; }\n"
        "  .title { font-size: 16px; font-weight: bold; }\n"
        "  .stats { font-size: 12px; fill: #555; }\n"
        "  .tick { font-size: 11px; }\n"
        "  .axis-label { font-size: 13px; }\n"
        "  .y-axis-label { font-size: 12px; }\n"
        "</style>"
    )
    lines.append('<rect width="100%" height="100%" fill="#fcfcfc"/>')

    lines.append(
        f'<text x="{svg_w // 2}" y="24" text-anchor="middle" class="title">'
        f"{chart_title} - All {n_subjects} Subjects "
        f"({len(all_vals)} bins, {zero_pct:.1f}% 0.000000)</text>"
    )
    lines.append(
        f'<text x="{svg_w // 2}" y="44" text-anchor="middle" class="stats">'
        f"mean={mn:.6f}  median={md:.6f}  IQR=[{q1:.6f}, {q3:.6f}]</text>"
    )

    lines.append(
        f'<line x1="{cx}" y1="{cy}" x2="{cx}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<line x1="{cx}" y1="{cy + chart_h}" '
        f'x2="{cx + chart_w}" y2="{cy + chart_h}" '
        f'stroke="#333" stroke-width="1"/>'
    )
    lines.append(
        f'<text x="18" y="{cy + chart_h // 2}" '
        f'text-anchor="middle" class="y-axis-label" '
        f'transform="rotate(-90 18,{cy + chart_h // 2})">'
        f"# bins (log)</text>"
    )

    bar_h_fn = _make_log_y_axis(lines, cx, cy, chart_h, max_count)

    n_nonzero_bins = num_buckets - 1
    bucket_width = (
        (max_val / n_nonzero_bins) if (max_val > 0 and n_nonzero_bins > 0) else 0.0
    )

    for b_idx in range(num_buckets):
        c = counts[b_idx]
        if c == 0:
            continue
        bh = bar_h_fn(c)
        bx = cx + b_idx * bar_w
        by = cy + chart_h - bh
        if b_idx == 0:
            tip = f"[0, 0]: {c}"
        else:
            rs = (b_idx - 1) * bucket_width
            re = b_idx * bucket_width
            tip = f"({rs:.6f}, {re:.6f}]: {c}"
        lines.append(
            f'<rect x="{bx + 0.5:.1f}" y="{by:.1f}" '
            f'width="{bar_w - 1:.1f}" height="{bh:.1f}" '
            f'fill="{bar_color}" fill-opacity="0.85">'
            f"<title>{tip}</title></rect>"
        )

    n_x_ticks = 5
    x_max = max_val if max_val > 0 else 1.0
    for t in range(n_x_ticks + 1):
        val = x_max * t / n_x_ticks
        x_pos = cx + chart_w * t / n_x_ticks
        lines.append(
            f'<text x="{x_pos:.1f}" y="{cy + chart_h + 16}" '
            f'text-anchor="middle" class="tick">{val:.6f}</text>'
        )
    lines.append(
        f'<text x="{cx + chart_w // 2}" y="{cy + chart_h + 38}" '
        f'text-anchor="middle" class="axis-label">{x_label}</text>'
    )
    lines.append("</svg>")

    with open(path, "w", encoding="utf-8") as fout:
        fout.write("\n".join(lines))


# ============================================================================
# Main
# ============================================================================


def main():
    print("daniel_path exists:", Path(daniel_path).exists())
    print("mira_path exists:  ", Path(mira_path).exists())
    print("deep_path exists:  ", Path(deep_path).exists())

    all_subjects, daniel_files, mira_files, deep_files = find_matching_files()

    # We need all three sources to compute all three pairs for a subject.
    common = [
        sid
        for sid in all_subjects
        if sid in daniel_files and sid in mira_files and sid in deep_files
    ]

    print(f"\nSubjects with all three sources: {len(common)} / {len(all_subjects)}")
    if not common:
        print("Nothing to compare.")
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # per_pair[pair_key] = list of (subject_id, [1-PCC, ...])
    per_pair = {
        "daniel_vs_mira": [],
        "daniel_vs_deep": [],
        "mira_vs_deep": [],
    }

    for sid in common:
        print(f"  {sid}...", end="", flush=True)
        try:
            daniel_list = read_daniel_template_mat(str(daniel_files[sid]))
            mira_list = read_mira_template_bin(str(mira_files[sid]))
            deep_data = read_deep_template_mat(str(deep_files[sid]))
        except Exception as e:
            print(f" FAILED ({e})")
            continue

        deep_list = deep_data.get("templates", [])
        deep_sr = deep_data.get("sample_rate", DEEP_SR_DEFAULT)

        dm, dp, mp = score_subject(daniel_list, mira_list, deep_list, deep_sr)
        per_pair["daniel_vs_mira"].append((sid, dm))
        per_pair["daniel_vs_deep"].append((sid, dp))
        per_pair["mira_vs_deep"].append((sid, mp))

        print(
            f"  Daniel-Mira n={len(dm)}, Daniel-Deep n={len(dp)}, Mira-Deep n={len(mp)}"
        )

    pretty = {
        "daniel_vs_mira": "Daniel vs Mira",
        "daniel_vs_deep": "Daniel vs Deep",
        "mira_vs_deep": "Mira vs Deep",
    }
    # Different bar colors so the three pairs are visually distinct when
    # you flip between SVGs.
    palette = {
        "daniel_vs_mira": ("#1a5fa8", "#2ecc71", "#3498db"),
        "daniel_vs_deep": ("#c0392b", "#9b59b6", "#e74c3c"),
        "mira_vs_deep": ("#16a085", "#f39c12", "#1abc9c"),
    }

    for key, vals in per_pair.items():
        if not vals:
            continue
        per_bar, ref_bar, pooled_bar = palette[key]
        title = f"1 - PCC between {pretty[key]} ECG Templates"

        per_path = OUTPUT_DIR / f"{key}_per_subject_histograms.svg"
        pooled_path = OUTPUT_DIR / f"{key}_pooled_histogram.svg"

        write_per_subject_pcc_grid(
            per_path,
            vals,
            chart_title=title,
            bar_color=per_bar,
            ref_color=ref_bar,
        )
        write_pooled_pcc_histogram(
            pooled_path,
            vals,
            chart_title=title,
            bar_color=pooled_bar,
        )

        print(f"  Wrote {per_path}")
        print(f"  Wrote {pooled_path}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    main()
