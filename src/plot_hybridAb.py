import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hybrid_in', type=str,
                        help='Path of the comet+hypedsearch output folder', required=True)
    parser.add_argument('--score', type=str,
                        help='Hybrid score (column) as the color bar', required=True)
    parser.add_argument('--outfile', type=str,
                        help='Name of the output file', required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    ab_hybrid_df = pd.read_csv(args.hybrid_in, index_col=0)
    ab_hybrid_df.sort_values(by=[args.score], ascending=True, inplace=True)
    ab_hybrid_df = ab_hybrid_df.reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(10,8))
    m = ax.scatter(x=ab_hybrid_df['left_abundance'],y=ab_hybrid_df['right_abundance'],c=ab_hybrid_df[args.score], cmap="plasma_r")
    ax.set_xlabel('Left Parent Abundance')
    ax.set_ylabel('Right Parent Abundance')
    n = ab_hybrid_df['name'].to_numpy()
    x = ab_hybrid_df['left_abundance'].to_numpy()
    y = ab_hybrid_df['right_abundance'].to_numpy()
    l = ab_hybrid_df[args.score].to_numpy()
    for i, txt in enumerate(n):
        if x[i] > 100 and y[i] > 100:
            ax.annotate(txt, (x[i]-len(txt)*1.25, y[i]+1), fontsize = 7)
        elif x[i] > 25 and y[i] > 25:
            ax.annotate(txt, (x[i]+1, y[i]+1), fontsize = 7)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    fig.colorbar(m, cax=cbar_ax,label=args.score)
    plt.savefig(args.outfile, dpi=fig.dpi, bbox_inches='tight')


if __name__ == '__main__':
    main()

