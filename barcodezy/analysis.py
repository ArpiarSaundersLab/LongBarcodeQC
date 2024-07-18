import pandas as pd
import parasail
import matplotlib.pyplot as plt

def read_length_hist(outpath: str, csv_name: str, alignment_counts: dict):
    csv = pd.read_csv(csv_name)
    experiment_name = os.path.basename()

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