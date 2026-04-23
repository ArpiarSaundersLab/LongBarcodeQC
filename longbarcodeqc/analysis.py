import os
import math
from io import StringIO
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from jinja2 import Environment, PackageLoader

def z_score_barcode_calling(reads_df: pd.DataFrame, z_thresh: float | None) -> Tuple[pd.DataFrame, str]:
    """Compute per-read barcode z-scores and top calls.

    Returns (summary_df, svg_plot_str).
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
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    axes[0].hist(z_scores[target_mask].values.flatten(), bins=35)
    axes[0].set_xlabel('z-score')
    axes[0].set_title('z-score distribution (linear)')
    axes[0].axvline(z_thresh, color='red', linestyle='--',
                    label='z-score threshold\nfor calling barcodes')
    axes[0].legend()

    axes[1].hist(z_scores[target_mask].values.flatten(), bins=35)
    axes[1].set_xlabel('z-score')
    axes[1].set_title('z-score distribution (log scale)')
    axes[1].set_yscale('log')
    axes[1].axvline(z_thresh, color='red', linestyle='--',
                    label='z-score threshold\nfor calling barcodes')
    axes[1].legend()

    plt.tight_layout()
    # save encoded plot
    svg_io = StringIO()
    plt.savefig(svg_io, format='svg', bbox_inches='tight')
    svg_content = svg_io.getvalue()
    svg_io.close()
    plt.close(fig)

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
    return summary_df, svg_content


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


def read_length_hist(
    summary_df: pd.DataFrame,
    hue_group_col: str,
    max_len: int,
    title: str,
    expected_insertions: int | None,
) -> str:
    """Generate a read length histogram and return it as an SVG string."""
    df_nofailed = summary_df[
        ~summary_df.read_type.str.endswith('_failed_anchor')
        ].copy()

    if hue_group_col == 'insertions':
    # add insertion group column (categorizes into max insertions vs < max insertions)
        num_pos = len([col for col in df_nofailed.columns if col.endswith('_1st_z_score')])
        if expected_insertions is not None:
            num_pos = expected_insertions
        # exclude ecoli reads from insertion histogram
        df_nofailed = df_nofailed[df_nofailed.read_type != 'ecoli']
        df_nofailed['ins_group'] = df_nofailed['insertions'] == num_pos
        df_nofailed['ins_group'] = df_nofailed.ins_group.replace({True: f'{num_pos}',
                                                                  False: f'0-{num_pos-1}'})
        if expected_insertions == 1:
            df_nofailed['ins_group'] = df_nofailed.ins_group.replace({f'0-{num_pos-1}': '0'})
        hue_group_col = 'ins_group'
    elif 'RE_' in hue_group_col:
        df_nofailed = df_nofailed[df_nofailed.read_type != 'ecoli']

    fig = plt.figure(figsize=(12, 7))
    sns.histplot(
        data=df_nofailed[df_nofailed.seq_len < max_len],
        x='seq_len',
        hue=hue_group_col,
        binwidth=_binwidth(max_len),
        alpha=0.5,
        palette='Set2',
        element='step',
    )
    sns.despine()
    plt.grid(visible=True, which='major', linestyle='--', alpha=0.4)
    plt.xlim(0, max_len)
    plt.xlabel('Read length (bp)')
    plt.title(title)
    legend = plt.gca().get_legend()
    if legend is not None:
        legend.set_title('')
    # save encoded plot
    svg_io = StringIO()
    plt.savefig(svg_io, format='svg', bbox_inches='tight')
    svg_content = svg_io.getvalue()
    svg_io.close()
    plt.close(fig)

    return svg_content

def MCS_length_hist(
    summary_df: pd.DataFrame,
    hue_group_col: str,
    title: str,
    expected_insertions: int | None,
) -> str:
    """Generate an MCS cassette length histogram and return it as an SVG string."""
    df_nofailed = summary_df[
        ~summary_df.read_type.str.endswith('_failed_anchor') &
        (summary_df.read_type != 'ecoli')
        ].copy()
    
    if hue_group_col == 'insertions':
    # add insertion group column (categorizes into max insertions vs < max insertions)
        num_pos = len([col for col in df_nofailed.columns if col.endswith('_1st_z_score')])
        if expected_insertions is not None:
            num_pos = expected_insertions
        # exclude ecoli reads from insertion histogram
        df_nofailed['ins_group'] = df_nofailed['insertions'] == num_pos
        df_nofailed['ins_group'] = df_nofailed.ins_group.replace({True: f'{num_pos}',
                                                                  False: f'0-{num_pos-1}'})
        if expected_insertions == 1:
            df_nofailed['ins_group'] = df_nofailed.ins_group.replace({f'0-{num_pos-1}': '0'})
        hue_group_col = 'ins_group'

    mcs_max = 2.5 * np.mean(df_nofailed.MCS_len)
    fig = plt.figure(figsize=(12, 7))
    sns.histplot(
        data=df_nofailed[df_nofailed.MCS_len < mcs_max],
        x='MCS_len',
        hue=hue_group_col,
        binwidth=_binwidth(mcs_max),
        alpha=0.5,
        palette='Set2',
        element='step',
    )
    sns.despine()
    plt.grid(visible=True, which='major', linestyle='--', alpha=0.4)
    plt.xlim(0, mcs_max)
    plt.xlabel('MCS cassette length (bp)')
    plt.title(title)
    legend = plt.gca().get_legend()
    if legend is not None:
        legend.set_title('')
    # save encoded plot
    svg_io = StringIO()
    plt.savefig(svg_io, format='svg', bbox_inches='tight')
    svg_content = svg_io.getvalue()
    svg_io.close()
    plt.close(fig)

    return svg_content

def report_gen(
    outpath: str,
    reads_df: pd.DataFrame,
    alignment_counts: Dict[str, int],
    ref_len: int,
    z_thresh: float | None,
    expected_insertions: int | None,
) -> pd.DataFrame:
    """Generate plots and HTML report, returning the summary DataFrame."""
    experiment_name = os.path.basename(outpath)

    ## alignment summary counts
    total_reads = sum(alignment_counts.values())
    align_df = pd.DataFrame(alignment_counts, index=['Count']).transpose()
    align_df = align_df.reset_index().rename(columns={'index': 'Alignment type'})
    percents = (align_df['Count'] / align_df['Count'].sum()) * 100
    percents = [f'{x:.2f}%' for x in percents]
    align_df['%'] = percents

    ## scoring table processing
    reads = reads_df.set_index('seq_id')
    # summary table output
    summary, z_plot = z_score_barcode_calling(reads, z_thresh)

    # get upper x lim for read length histograms
    max_len = math.ceil(ref_len / 1000) * 1000
    max_len = max_len + math.ceil(max_len*0.1 / 1000) * 1000 # add ~10% buffer to max_len

    read_type_hist = read_length_hist(summary, 'read_type', max_len, 
                                      'Read type', expected_insertions)
    insertion_hist = read_length_hist(summary, 'insertions', max_len, 
                                      'Number of valid position insertions', 
                                      expected_insertions)
    MCS_insertion_hist = MCS_length_hist(summary, 'insertions', 
                                        'Insertions by cassette length', expected_insertions)

    # plots for each restriction site
    rs_plots = []
    rs_mcs_plots = []
    for rs in summary.loc[:, summary.columns.str.contains('RE_')].columns:
        rs_plots.append(read_length_hist(summary, rs, max_len,
                                         f'{rs[3:]} site detection', expected_insertions))
        rs_mcs_plots.append(MCS_length_hist(summary, rs,
                                            f'{rs[3:]} site detection (MCS length)', expected_insertions))
    
    # restriction site % table
    enz_names = summary.loc[:, summary.columns.str.contains('RE_')].columns
    df_no_ecoli_nofailed = summary[
        ~summary.read_type.str.endswith('_failed_anchor') &
        (summary.read_type != 'ecoli')
        ].copy()
    # avoid division by zero if no plasmid reads in expected range
    if df_no_ecoli_nofailed.shape[0] == 0:
        print('Warning: No non-bacterial reads found.')
        enz_all = ['N/A' for _ in enz_names]
    else:
        enz_all = [f'{100*sum(df_no_ecoli_nofailed[x])/df_no_ecoli_nofailed.shape[0]:.2f}%' for x in enz_names]
    # In default case (AP-Amp/AP-Kan), all non-ecoli/non-failed reads are the target.
    # In custom plasmid case, filter to 'plasmid' read_type only.
    if 'plasmid' in df_no_ecoli_nofailed.read_type.values:
        target_mask = df_no_ecoli_nofailed.read_type == 'plasmid'
    else:
        target_mask = pd.Series(True, index=df_no_ecoli_nofailed.index)
    df_plasmid_ranged = df_no_ecoli_nofailed[
        target_mask &
        (df_no_ecoli_nofailed.seq_len > int(ref_len*0.85)) &
        (df_no_ecoli_nofailed.seq_len < int(ref_len*1.15))
        ]
    # avoid division by zero if no plasmid reads in expected range
    if df_plasmid_ranged.shape[0] == 0:
        print('Warning: No plasmid reads found in expected length range.')
        enz_specific = ['N/A' for _ in enz_names]
    else:
        enz_specific = [f'{100*sum(df_plasmid_ranged[x])/df_plasmid_ranged.shape[0]:.2f}%' for x in enz_names] 

    enz_df = pd.DataFrame({'Restriction site': [x[3:] for x in enz_names],
                           f'% in MCS sequence for plasmid of interest': enz_specific, 
                           f'% in full length non-bacterial reads': enz_all})


    # generate html 
    env = Environment(loader=PackageLoader('longbarcodeqc', 'template'))
    template = env.get_template('template.html')

    context = {
     'title': experiment_name,
     'subtitle': 'Library summary report',
     'headers': align_df.columns.tolist(),
     'data': align_df.values.tolist(),
     'total_reads': total_reads,
     'read_type': read_type_hist,
     'insertions': insertion_hist,
     'rs_plots': [{'read': r, 'mcs': m} for r, m in zip(rs_plots, rs_mcs_plots)],
     'rs_plot_names': enz_df['Restriction site'].tolist(),
     'headers_enz': enz_df.columns.tolist(),
     'data_enz': enz_df.values.tolist(),
     'z_plot': z_plot,
     'MCS_insertions': MCS_insertion_hist
    }

    html_report = template.render(context)

    with open(f'{outpath}/{experiment_name}_summary_report.html', 'w') as fh:
        fh.write(html_report)
    
    return summary
