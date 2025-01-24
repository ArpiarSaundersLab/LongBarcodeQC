import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from jinja2 import Environment, PackageLoader, select_autoescape
import base64
from io import StringIO
import math

def z_score_barcode_calling(reads_df, z_thresh=5.5):
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
    z_scores = (reads_df[site_columns].
                sub(reads_df[site_columns].mean(axis=1), axis=0).
                div(reads_df[site_columns].std(axis=1), axis=0).
                fillna(0))

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
    return summary_df


def read_length_hist(summary_df, hue_group_col: str, max_len: int, title: str):
    df_nofailed = summary_df[
        (summary_df.read_type != 'plasmid_failed_anchor') &
        (summary_df.read_type != 'a31_failed_anchor')
        ].copy()
    
    if hue_group_col == 'insertions':
    # add insertion group column (categorizes into max insertions vs < max insertions)
        num_pos = len([col for col in df_nofailed.columns if col.endswith('_1st_z_score')])
        df_nofailed['ins_group'] = df_nofailed['insertions'] == num_pos
        df_nofailed['ins_group'] = df_nofailed.ins_group.replace({True: f'{num_pos}',
                                                                  False: f'0-{num_pos-1}'})
        hue_group_col = 'ins_group'

    plt.figure(figsize=(12, 7))
    sns.histplot(data=df_nofailed[df_nofailed.seq_len < max_len], x='seq_len', 
                 hue=hue_group_col, bins=60, alpha=0.5, palette='Set2', element='step')
    sns.despine()
    plt.grid(visible=True, which='major', linestyle='--', alpha=0.4)
    plt.xlim(0, max_len)
    plt.xlabel('Read length (bp)')
    plt.title(title)
    plt.gca().get_legend().set_title('')
    # save encoded plot
    svg_io = StringIO()
    plt.savefig(svg_io, format='svg', bbox_inches='tight')
    svg_content = svg_io.getvalue()
    svg_io.close()

    return svg_content

def report_gen(outpath: str, reads_df, alignment_counts: str, ref_len: int):
    experiment_name = os.path.basename(outpath)

    ## alignment summary counts
    total_reads = sum(alignment_counts.values())
    #alignment_counts['Total'] = sum(alignment_counts.values())
    align_df = pd.DataFrame(alignment_counts, index=['Count']).transpose()
    align_df = align_df.reset_index().rename(columns={'index': 'Alignment type'})
    percents = (align_df['Count'] / align_df['Count'].sum()) * 100
    percents = [f'{x:.2f}%' for x in percents]
    align_df['%'] = percents

    ## scoring table processing
    reads = reads_df.set_index('seq_id')
    # summary table output
    summary = z_score_barcode_calling(reads)

    # get upper x lim for read length histograms
    max_len = math.ceil(ref_len / 1000) * 1000
    max_len = max_len + math.ceil(max_len*0.1 / 1000) * 1000 # add ~10% buffer to max_len

    read_type_hist = read_length_hist(summary, 'read_type', max_len, 
                                      'Read type')
    insertion_hist = read_length_hist(summary, 'insertions', max_len, 
                                      'Number of valid position insertions')

    # generate html 
    env = Environment(loader=PackageLoader('barcodezy', 'template'))
    template = env.get_template('template.html')

    context = {
     'title': experiment_name,
     'subtitle': 'Library summary report',
     'headers': align_df.columns.tolist(),
     'data': align_df.values.tolist(),
     'total_reads': total_reads,
     'read_type': read_type_hist,
     'insertions': insertion_hist
    }

    html_report = template.render(context)

    with open(f'{outpath}/{experiment_name}_summary_report.html', 'w') as fh:
        fh.write(html_report)
    
    return summary