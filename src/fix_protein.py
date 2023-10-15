import sys
import os
import argparse
import numpy as np
import pandas as pd
sys.path.insert(0, 'src')
import utils as ut


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comet_in', type=str,
                        help='Path of the comet+hypedsearch output folder', required=True)
    parser.add_argument('--hypedsearch_in', type=str,
                        help='Path of the hypedsearch output folder', required=True)
    parser.add_argument('--out', type=str,
                        help='Path of the output file', required=True)
    args = parser.parse_args()
    return args


def fix_hybrid(df_hybrid_cut):
    '''
    Remove redundancy in hybrid-parent data (Hypedsearch_outputs)
    '''
    df_hybrid = df_hybrid_cut.copy()
    df_hybrid_cut['right kmer'] = df_hybrid['right kmer'].str.replace('/', ',', regex=False).str.split(",")
    df_hybrid_cut = df_hybrid_cut.explode('right kmer')
    df_hybrid_cut['left kmer'] = df_hybrid['left kmer'].str.replace('/', ',', regex=False).str.split(",")
    df_hybrid_cut = df_hybrid_cut.explode('left kmer')
    df_hybrid_cut = df_hybrid_cut.drop_duplicates(subset=['sequence','left kmer','right kmer'])
    df_hybrid_cut_r = df_hybrid_cut.groupby(['sequence'])['right kmer'].apply(lambda x: '/'.join(set(x))).reset_index()
    df_hybrid_cut_l = df_hybrid_cut.groupby(['sequence'])['left kmer'].apply(lambda x: '/'.join(set(x))).reset_index()
    df_hybrid_cut_fixed = df_hybrid_cut_l.merge(df_hybrid_cut_r,how='inner')
    return df_hybrid_cut_fixed


def fix_run(df, df_hybrid_cut_fixed):
    '''
    Fix Comet+Hypedsearch result matching hybrid-parent data
    '''
    ind_list = df.index
    protein = df['protein'].str.split(',', expand=True)
    protein = protein.astype(str)
    seq = []
    for i in range (len(protein)):
        find = protein.iloc[i].str.find('Hypedsearch_Hybrid_')
        for j in np.where(find!=-1)[0]:
            seq.append(protein.iloc[i][j][find[j]+19:])
    sequances = pd.DataFrame({'sequence':seq})
    sequances_merged = sequances.merge(df_hybrid_cut_fixed, on ='sequence',how='left',indicator=True)[['sequence','left kmer', 'right kmer','_merge']]
    protein_fixed = protein.copy()
    k = 0
    for i in range (len(protein_fixed)):
        find = protein_fixed.iloc[i].str.find('Hypedsearch_Hybrid_')
        for j in np.where(find!=-1)[0]:
            protein_fixed.iloc[i][j] = sequances_merged['left kmer'].iloc[k]+'-'+sequances_merged['right kmer'].iloc[k]+protein_fixed.iloc[i][j][find.iloc[j]-1:]
            k = k +1
    protein_fixed_pd = pd.DataFrame(protein_fixed)
    protein_fixed_pd.replace(to_replace='None', value=np.nan, inplace=True)
    protein_fixed_pd.index = ind_list
    df['protein_fix'] = protein_fixed_pd[protein_fixed_pd.columns[:]].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1)
    return df


def main():
    args = get_args()
    df_hybrid = ut.getData(args.hypedsearch_in, skiprow1=False)
    df_hybrid_cut_fixed = fix_hybrid(df_hybrid.copy())
    df_run = ut.getData(args.comet_in)
    df_run_nodecoy = df_run[~df_run.protein.str.contains("DECOY")]
    df_run_nodecoy_fixed = fix_run(df_run_nodecoy.copy(), df_hybrid_cut_fixed.copy())
    df_run_nodecoy_fixed.to_csv(os.path.join(args.out,"fixed_comethypedsearch_outputs.csv"))


if __name__ == '__main__':
    main()
