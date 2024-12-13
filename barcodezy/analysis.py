import os
import pandas as pd
import matplotlib.pyplot as plt
from jinja2 import Environment, PackageLoader, select_autoescape
import base64
from io import BytesIO

def read_length_hist(outpath: str):
    experiment_name = os.path.basename(outpath)
    csv = pd.read_csv(f'{outpath}/{experiment_name}.csv')

    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(7,4))
    ax.hist(csv.seq_len, bins=30, edgecolor='black')
    ax.set_xlabel('Read length (bp)')
    ax.set_ylabel('count')
    ax.set_title('Read length distribution')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.grid(alpha=0.4)
    fig.savefig(f'{outpath}/{experiment_name}_read_length_distribution.png')

def report_gen(outpath: str, reads_df, alignment_counts: str):
    experiment_name = os.path.basename(outpath)
    df = reads_df

    # alignment summary counts
    alignment_counts['Total'] = sum(alignment_counts.values())
    reads = pd.DataFrame(alignment_counts, index=['Count']).transpose()
    reads = reads.reset_index().rename(columns={'index': 'Alignment type'})

    # read length histogram
    plt.style.use('seaborn-v0_8-pastel')
    fig, ax = plt.subplots(figsize=(9,4))
    ax.hist(df.seq_len, bins=40, edgecolor='black')
    ax.set_xlabel('Read length (bp)')
    ax.set_ylabel('count')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(left=0)
    ax.grid(alpha=0.4)

    buf = BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()

    env = Environment(loader=PackageLoader('barcodezy', 'template'))
    template = env.get_template('template.html')

    context = {
     'title': experiment_name,
     'subtitle': 'Library summary report',
     'headers': reads.columns.tolist(),
     'data': reads.head().values.tolist(),
     'img_base64': img_base64
    }

    html_report = template.render(context)

    with open(f'{outpath}/{experiment_name}_summary_report.html', 'w') as fh:
        fh.write(html_report)