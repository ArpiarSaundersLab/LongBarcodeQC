import os
import math
from datetime import datetime
from io import StringIO
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from jinja2 import Environment, PackageLoader

# report plot theme: quiet chrome so the data carries the color
_INK = '#0b0b0b'
_INK_2 = '#52514e'
_MUTED = '#898781'
_GRID = '#e1e0d9'
_AXIS = '#c3c2b7'
_SURFACE = '#ffffff'
_NEUTRAL = '#b5b3ab'
# categorical slots, assigned in this order (validated for color-vision deficiency)
_SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']

_THEME = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'text.color': _INK_2,
    'axes.edgecolor': _AXIS,
    'axes.labelcolor': _INK_2,
    'axes.labelsize': 10,
    'axes.linewidth': 1,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': False,
    'axes.grid': True,
    'axes.grid.axis': 'y',
    'axes.axisbelow': True,
    'axes.facecolor': _SURFACE,
    'figure.facecolor': _SURFACE,
    'grid.color': _GRID,
    'grid.linewidth': 1,
    'grid.linestyle': '-',
    'xtick.color': _MUTED,
    'ytick.color': _MUTED,
    'xtick.labelcolor': _MUTED,
    'ytick.labelcolor': _MUTED,
    'ytick.left': False,
    'legend.frameon': False,
    'legend.fontsize': 10,
    'svg.fonttype': 'path',
}

_thousands = FuncFormatter(lambda x, _: f'{x:,.0f}')


def _svg(fig) -> str:
    svg_io = StringIO()
    fig.savefig(svg_io, format='svg', bbox_inches='tight')
    plt.close(fig)
    return svg_io.getvalue()


def _top_legend(ax, ncol: int) -> None:
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=ncol, handlelength=1.2,
              columnspacing=1.6, borderaxespad=0.2, labelcolor=_INK_2)


def _target_type(reads_df: pd.DataFrame) -> str:
    """Read type of the target plasmid: 'plasmid' for a custom/SBARRO plasmid, else 'AP-Amp'."""
    return 'plasmid' if reads_df.read_type.str.startswith('plasmid').any() else 'AP-Amp'


def _read_type_style(read_types: list[str], target: str) -> Dict[str, tuple[str, str]]:
    """Map each read type to (legend label, color). Color follows the read type, not its rank."""
    def label(read_type: str) -> str:
        base = read_type.replace('_failed_anchor', '')
        name = {'plasmid': 'Target plasmid', 'ecoli': 'E. coli'}.get(base, base)
        if base == target and target != 'plasmid':
            name = f'{base} (target)'
        return f'{name}, MCS not found' if read_type.endswith('_failed_anchor') else name

    order = [target, f'{target}_failed_anchor']
    order += sorted(t for t in read_types if t not in order and t != 'ecoli')
    styles = {t: (label(t), _SERIES[i % len(_SERIES)]) for i, t in enumerate(order)}
    styles['ecoli'] = (label('ecoli'), _NEUTRAL)
    return {t: styles[t] for t in order + ['ecoli'] if t in read_types}


def z_score_barcode_calling(reads_df: pd.DataFrame,
                            z_thresh: float | None) -> Tuple[pd.DataFrame, str, float]:
    """Compute per-read barcode z-scores and top calls.

    Returns (summary_df, svg_plot_str, z-score threshold used).
    """
    error_msg = ('Could not parse barcode names. Naming scheme must be "siteX_posY_Z."\n'
                 'X: site number, Y: = position, Z:  barcode #\n'
                 'e.g. site2_posB_127')

    site_columns = [col for col in reads_df if col.startswith('site')]

    # error if no columns start with "site":
    if len(site_columns) < 1: raise ValueError(error_msg)
    # error if site names don't follow site_pos_bc scheme:
    for bc_name in site_columns:
        split_name = bc_name.split('_')
        if len(split_name) != 3: raise ValueError(error_msg)

    # obtain unique site_pos groups:
    site_pos_groups = set([f'{x.split("_")[0]}_{x.split("_")[1]}' for x in site_columns])
    site_pos_groups = sorted(list(site_pos_groups))

    # z score matrix of BC alignment scores
    # Convert the matrix to a NumPy array and ignore NaNs explicitly
    global_mean = np.nanmean(reads_df[site_columns].values)
    global_std = np.nanstd(reads_df[site_columns].values)

    # Compute z-scores using NumPy (element-wise)
    z_scores = (reads_df[site_columns] - global_mean) / global_std

    # Fill NaNs (if any) with 0
    z_scores = z_scores.fillna(0)

    # Target reads mask: 'plasmid' for custom plasmid, or all non-ecoli/non-failed for default AP case
    if 'plasmid' in reads_df.read_type.values:
        target_mask = reads_df.read_type == 'plasmid'
    else:
        target_mask = ~reads_df.read_type.str.endswith('_failed_anchor') & (reads_df.read_type != 'ecoli')

    # Use minimum z-score as threshold for calling insertions if not provided
    if z_thresh is None:
        z_thresh = abs(min(z_scores[target_mask].values.flatten()))

    # z-score distribution plots (linear and log scaled)
    values = z_scores[target_mask].values.flatten()
    with plt.rc_context(_THEME):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
        for ax, scale in zip(axes, ('linear', 'log')):
            ax.hist(values, bins=35, color=_SERIES[0], edgecolor=_SURFACE, linewidth=0.8)
            ax.set_yscale(scale)
            if scale == 'linear':
                ax.yaxis.set_major_formatter(_thousands)
            ax.axvline(z_thresh, color=_INK, linestyle=(0, (4, 3)), linewidth=1.2)
            # label on whichever side of the line has room
            right_side = z_thresh < values.min() + 0.7 * np.ptp(values)
            ax.annotate(f'threshold {z_thresh:.2f}', xy=(z_thresh, 1), xycoords=('data', 'axes fraction'),
                        xytext=(6 if right_side else -6, -4), textcoords='offset points', va='top',
                        ha='left' if right_side else 'right', color=_INK, fontsize=9,
                        bbox=dict(boxstyle='square,pad=0.15', fc=_SURFACE, ec='none'))
            ax.set_xlabel('z-score')
            ax.set_title(f'{scale.capitalize()} scale', loc='left', fontsize=10, color=_INK_2)
        axes[0].set_ylabel('Barcode alignments')
        fig.tight_layout()
        svg_content = _svg(fig)

    # obtain top alignment for each site_pos group (e.g. site1_posA)
    top_scores_per_group = {}
    for group in site_pos_groups:
        top_scores_per_group[f'{group}_1st_name'] = (z_scores.transpose().
                                                     filter(like=group, axis=0).
                                                     idxmax().str.split('_', expand=True)[2])
        top_scores_per_group[f'{group}_1st_z_score'] = (z_scores.transpose().
                                                        filter(like=group, axis=0).max())
    barcode_calls = pd.DataFrame.from_dict(top_scores_per_group)

    insertions = (barcode_calls[[f'{x}_1st_z_score' for x in site_pos_groups]] > z_thresh).sum(axis=1)
    barcode_calls['insertions'] = insertions

    # add z_score barcode calls to main df
    summary_df = pd.concat([reads_df.drop(site_columns, axis=1), barcode_calls], axis=1)
    return summary_df, svg_content, float(z_thresh)


def _binwidth(max_val: float) -> int:
    if max_val < 100:
        return 2
    elif max_val < 1000:
        return 5
    elif max_val < 5000:
        return 20
    elif max_val < 10000:
        return 40
    return 80


def _num_positions(df: pd.DataFrame, expected_insertions: int | None) -> int:
    if expected_insertions is not None:
        return expected_insertions
    return len([col for col in df.columns if col.endswith('_1st_z_score')])


def _group_styles(df: pd.DataFrame, hue_group_col: str, expected_insertions: int | None,
                  target: str) -> tuple[pd.Series, Dict[str, tuple[str, str]]]:
    """Return (group per read, {group: (legend label, color)}) for a histogram hue."""
    if hue_group_col == 'read_type':
        return df.read_type, _read_type_style(df.read_type.unique().tolist(), target)
    if hue_group_col == 'insertions':
        num_pos = _num_positions(df, expected_insertions)
        full = df['insertions'] == num_pos
        partial = '0' if num_pos == 1 else f'0–{num_pos - 1}'
        groups = full.map({True: 'full', False: 'partial'})
        return groups, {'full': (f'{num_pos} positions called', _SERIES[0]),
                        'partial': (f'{partial} positions called', _SERIES[1])}
    # restriction site columns
    groups = df[hue_group_col].astype(bool).map({True: 'present', False: 'absent'})
    return groups, {'present': ('Site present', _SERIES[0]), 'absent': ('Site absent', _NEUTRAL)}


def _length_hist(values: pd.Series, groups: pd.Series, styles: Dict[str, tuple[str, str]],
                 max_len: float, xlabel: str) -> str:
    """Overlaid step histograms (one per group) with a light area wash."""
    edges = np.arange(0, max_len + _binwidth(max_len), _binwidth(max_len))
    with plt.rc_context(_THEME):
        fig, ax = plt.subplots(figsize=(12, 4.2))
        for group, (label, color) in styles.items():
            counts, _ = np.histogram(values[groups == group], bins=edges)
            if counts.sum() == 0:
                continue
            ax.stairs(counts, edges, fill=True, color=color, alpha=0.12, linewidth=0)
            ax.stairs(counts, edges, color=color, linewidth=1.5, label=label)
        ax.set_xlim(0, max_len)
        ax.set_ylim(bottom=0)
        ax.xaxis.set_major_formatter(_thousands)
        ax.yaxis.set_major_formatter(_thousands)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Reads')
        _top_legend(ax, ncol=len(styles))
        return _svg(fig)


def read_length_hist(
    summary_df: pd.DataFrame,
    hue_group_col: str,
    max_len: int,
    expected_insertions: int | None,
) -> str:
    """Generate a read length histogram and return it as an SVG string."""
    target = _target_type(summary_df)
    df = summary_df.copy()
    if hue_group_col != 'read_type':
        # barcode and restriction site groups are only defined for anchored plasmid reads
        df = df[~df.read_type.str.endswith('_failed_anchor') & (df.read_type != 'ecoli')]
    df = df[df.seq_len < max_len]
    groups, styles = _group_styles(df, hue_group_col, expected_insertions, target)
    return _length_hist(df.seq_len, groups, styles, max_len, 'Read length (bp)')


def MCS_length_hist(
    summary_df: pd.DataFrame,
    hue_group_col: str,
    expected_insertions: int | None,
) -> str:
    """Generate an MCS cassette length histogram and return it as an SVG string."""
    target = _target_type(summary_df)
    df = summary_df[
        ~summary_df.read_type.str.endswith('_failed_anchor') &
        (summary_df.read_type != 'ecoli')
        ].copy()
    mcs_max = 2.5 * np.mean(df.MCS_len)
    df = df[df.MCS_len < mcs_max]
    groups, styles = _group_styles(df, hue_group_col, expected_insertions, target)
    return _length_hist(df.MCS_len, groups, styles, mcs_max, 'MCS cassette length (bp)')


def backbone_map(backbone: dict, mcs_len: float) -> str:
    """Circular plasmid map of the backbone consensus with read depth around it.

    The MCS sits at the top (drawn to scale from the median MCS length), with the left
    flank before it and the right flank after it going clockwise; the backbone consensus
    runs clockwise from the right flank around to the left flank.
    """
    depth = backbone['depth']
    bb_len = len(depth)
    left_len, right_len = backbone['left_flank_len'], backbone['right_flank_len']
    total = bb_len + mcs_len

    def theta(pos):
        # position 0 = start of the backbone (right flank), just clockwise of the MCS
        return np.pi / 2 - 2 * np.pi * (mcs_len / 2 + np.asarray(pos, dtype=float)) / total

    ring, band = 1.0, 0.32
    with plt.rc_context(_THEME):
        fig = plt.figure(figsize=(6.4, 6.4))
        ax = fig.add_subplot(projection='polar')
        ax.set_axis_off()
        ax.set_ylim(0, ring + band + 0.42)

        # read depth as a band outside the backbone ring; a rolling median removes the
        # single-base dips from read indels (mostly homopolymers), which are not coverage loss
        window = max(5, bb_len // 300)
        smooth = pd.Series(depth).rolling(window, center=True, min_periods=1).median().to_numpy()
        step = max(1, bb_len // 1500)
        pos = np.arange(0, bb_len, step)
        r = ring + 0.03 + band * smooth[pos] / max(smooth.max(), 1)
        ax.fill_between(theta(pos), ring + 0.03, r, color=_SERIES[0], alpha=0.14, linewidth=0)
        ax.plot(theta(pos), r, color=_SERIES[0], linewidth=1.2)

        # backbone ring, MCS arc and flank arcs
        ax.plot(theta(np.linspace(0, bb_len, 600)), np.full(600, ring), color=_AXIS, linewidth=3,
                solid_capstyle='butt')
        mcs = np.linspace(-mcs_len, 0, 50)
        ax.plot(theta(mcs), np.full(50, ring), color=_INK_2, linewidth=3, linestyle=(0, (1.5, 1.5)))
        flanks = [('Right flank', 0, right_len, _SERIES[1]),
                  ('Left flank', bb_len - left_len, bb_len, _SERIES[2])]
        for _, start, end, color in flanks:
            ax.plot(theta(np.linspace(start, end, 30)), np.full(30, ring), color=color, linewidth=7,
                    solid_capstyle='butt')

        # differences from the reference, as ticks just inside the ring
        variants = backbone['comparison']['variants'] if backbone.get('comparison') else []
        for v in variants:
            t = theta(v['consensus_pos'])
            ax.plot([t, t], [ring - 0.09, ring - 0.03], color=_INK, linewidth=1)

        # kb ticks inside the ring
        tick_step = next(s for s in (500, 1000, 2000, 5000, 10000, 20000) if bb_len / s <= 10)
        for kb in range(tick_step, bb_len - tick_step // 3, tick_step):
            t = theta(kb)
            ax.plot([t, t], [ring - 0.02, ring - 0.0], color=_AXIS, linewidth=1)
            ax.text(t, ring - 0.17, f'{kb / 1000:g} kb', ha='center', va='center',
                    color=_MUTED, fontsize=8)

        # labels for the MCS and flanks, pushed outward with short leader lines
        labels = [('MCS', -mcs_len / 2, _INK_2, f'~{mcs_len:,.0f} bp'),
                  ('Right flank', right_len / 2, _INK_2, f'{right_len} bp'),
                  ('Left flank', bb_len - left_len / 2, _INK_2, f'{left_len} bp')]
        spread = {'MCS': 0.0, 'Right flank': -0.32, 'Left flank': 0.32}
        for name, p, color, sub in labels:
            t = float(theta(p))
            tl = t + spread[name]
            ax.plot([t, tl], [ring + 0.04, ring + band + 0.16], color=_AXIS, linewidth=0.8)
            ha = 'center' if name == 'MCS' else ('left' if np.cos(tl) > 0 else 'right')
            ax.text(tl, ring + band + 0.22, f'{name}\n', ha=ha, va='center', color=_INK,
                    fontsize=10, fontweight='bold')
            ax.text(tl, ring + band + 0.22, f'\n{sub}', ha=ha, va='center', color=_MUTED, fontsize=9)

        # centre summary
        ax.text(0, 0, f'{bb_len:,} bp', ha='center', va='bottom', color=_INK, fontsize=20,
                fontweight='bold')
        ax.text(0, 0, f'\nbackbone consensus\n{backbone["reads_used"]:,} reads · '
                      f'mean depth {backbone["mean_depth"]:,.0f}×',
                ha='center', va='top', color=_INK_2, fontsize=9.5, linespacing=1.5)
        return _svg(fig)


def _headline_metrics(summary: pd.DataFrame, alignment_counts: Dict[str, int],
                      expected_insertions: int | None) -> dict:
    """Key numbers for the top of the report."""
    total = sum(alignment_counts.values())
    target = _target_type(summary)
    target_reads = summary[summary.read_type.isin([target, f'{target}_failed_anchor'])]
    anchored = target_reads[target_reads.read_type == target]
    num_pos = _num_positions(summary, expected_insertions)
    full = int((anchored.insertions == num_pos).sum())

    def pct(n, d):
        return 100 * n / d if d else 0.0

    return {
        'total_reads': total,
        'target_label': 'AP-Amp' if target == 'AP-Amp' else 'Target plasmid',
        'target_reads': len(target_reads),
        'target_pct': pct(len(target_reads), total),
        'anchored_reads': len(anchored),
        'anchored_pct': pct(len(anchored), len(target_reads)),
        'num_positions': num_pos,
        'full_reads': full,
        'median_positions': float(anchored.insertions.median()) if len(anchored) else float('nan'),
        'full_pct': pct(full, len(anchored)),
        'median_mcs_len': float(anchored.MCS_len.median()) if len(anchored) else float('nan'),
        'median_read_len': float(target_reads.seq_len.median()) if len(target_reads) else float('nan'),
    }


def report_gen(
    outpath: str,
    reads_df: pd.DataFrame,
    alignment_counts: Dict[str, int],
    z_thresh: float | None,
    expected_insertions: int | None,
    user_command: str,
    backbone: dict | None = None,
    version: str = '',
) -> pd.DataFrame:
    """Generate plots and HTML report, returning the summary DataFrame."""
    experiment_name = os.path.basename(outpath)

    ## alignment summary counts
    total_reads = sum(alignment_counts.values())
    align_rows = [(name, count, 100 * count / total_reads if total_reads else 0.0)
                  for name, count in alignment_counts.items()]

    ## scoring table processing
    reads = reads_df.set_index('seq_id')
    # summary table output
    summary, z_plot, z_thresh_used = z_score_barcode_calling(reads, z_thresh)

    # get upper x lim for read length histograms from actual read length distribution
    p99 = summary.seq_len.quantile(0.95)
    max_len = math.ceil(p99 / 1000) * 1000
    max_len = max_len + math.ceil(max_len * 0.1 / 1000) * 1000  # add ~10% buffer

    read_type_hist = read_length_hist(summary, 'read_type', max_len, expected_insertions)
    insertion_hist = read_length_hist(summary, 'insertions', max_len, expected_insertions)
    MCS_insertion_hist = MCS_length_hist(summary, 'insertions', expected_insertions)

    # restriction site % table and plots — based on reads with successful MCS anchoring only
    enz_names = summary.loc[:, summary.columns.str.contains('RE_')].columns
    df_mcs_anchored = summary[
        ~summary.read_type.str.endswith('_failed_anchor') &
        (summary.read_type != 'ecoli')
        ].copy()
    if df_mcs_anchored.shape[0] == 0:
        print('Warning: No reads with successful MCS anchoring found.')
    restriction_sites = []
    for rs in enz_names:
        pct = (100 * df_mcs_anchored[rs].astype(bool).sum() / df_mcs_anchored.shape[0]
               if df_mcs_anchored.shape[0] else None)
        restriction_sites.append({
            'name': rs[3:],
            'pct': pct,
            'read_plot': read_length_hist(summary, rs, max_len, expected_insertions),
            'mcs_plot': MCS_length_hist(summary, rs, expected_insertions),
        })

    metrics = _headline_metrics(summary, alignment_counts, expected_insertions)
    if backbone is not None:
        backbone = dict(backbone, map_plot=backbone_map(backbone, metrics['median_mcs_len']))

    # generate html
    env = Environment(loader=PackageLoader('longbarcodeqc', 'template'))
    template = env.get_template('template.html')

    context = {
        'title': experiment_name,
        'version': version,
        'generated': datetime.now().strftime('%Y-%m-%d %H:%M'),
        'metrics': metrics,
        'align_rows': align_rows,
        'read_type': read_type_hist,
        'insertions': insertion_hist,
        'MCS_insertions': MCS_insertion_hist,
        'z_plot': z_plot,
        'z_thresh': z_thresh_used,
        'z_thresh_user_set': z_thresh is not None,
        'restriction_sites': restriction_sites,
        'backbone': backbone,
        'user_command': user_command,
    }

    html_report = template.render(context)

    with open(f'{outpath}/{experiment_name}_summary_report.html', 'w') as fh:
        fh.write(html_report)

    return summary
