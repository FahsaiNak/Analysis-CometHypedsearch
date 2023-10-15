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


def countAb(df):
    protein_dict = dict()
    for cluster in df.protein_fix:
        cluster = ut.dropDup(cluster.split(','))
        proname_list = list()
        code = ['sp|', 'tr|', 'Hybrid'] #order sensitive
        for c in code:
            name_list = [name for name in cluster if c in name]
            for name in name_list:
                if name not in proname_list:
                    if c == 'Hybrid':
                        name = name.split('|')[0]
                    proname_list.append(name)
            if len(name_list) != 0:
                break
        for proname in proname_list:
            if proname not in protein_dict:
                protein_dict[proname] = 1
            else:
                protein_dict[proname] += 1
    return dict(sorted(protein_dict.items(), key=lambda item: item[1], reverse=True))


def createAbDf(df_run_nodecoy_fixed, top_protein_hybrid, top_protein_native, df_hybrid):
    ab_hybrid={'name':[],'hybrid_abundance':[],'left_abundance':[],'right_abundance':[],"max_ions_matched":[],"candidate_sequences":[],"all_sequences":[]}
    for rk, proname in enumerate(list(top_protein_hybrid.keys())):
        left = proname.split("-")[0]
        right = proname.split("-")[-1]
        df = df_run_nodecoy_fixed[(df_run_nodecoy_fixed.protein_fix.str.contains(left+"-"+right)) & (~df_run_nodecoy_fixed.protein_fix.str.contains("sp|tr"))]
        score = max(df.ions_matched)
        seq_list = list()
        lead_list = list()
        for ind in df.index:
            name = df.protein_fix[ind]
            for protein in ut.dropDup(name.split(",")):
                if protein.split("|Hypedsearch_Hybrid_")[0] == left+"-"+right:
                    seq = protein.split("|Hypedsearch_Hybrid_")[1]
                    seq_list.append(seq)
                    lead = ut.checkSeq(seq, 35, 'D', 'D', df_hybrid)
                if lead != "":
                    lead_list.append(lead)
        seqs = ",".join(ut.dropDup(seq_list))
        lead_seqs = ",".join(ut.dropDup(lead_list))
        if seqs != "":
            print(rk+1,"/",len(list(top_protein_hybrid.keys())), left+"-"+right, lead_seqs)
            ab_hybrid['name'].append(proname)
            ab_hybrid['hybrid_abundance'].append(top_protein_hybrid[proname])
            ab_hybrid['max_ions_matched'].append(score)
            ab_hybrid['candidate_sequences'].append(lead_seqs)
            ab_hybrid['all_sequences'].append(seqs)
            ab_hybrid['left_abundance'].append(ut.getNativeAb(left, top_protein_native))
            ab_hybrid['right_abundance'].append(ut.getNativeAb(right, top_protein_native))
            print("-------")
    return ab_hybrid


def main():
    args = get_args()
    df = pd.read_csv(args.comet_in, index_col=0)
    df_run_nodecoy_fixed = df[(df.num == 1) & (df.ions_matched >= 3)]
    df_hybrid = ut.getData(args.hypedsearch_in, skiprow1=False)
    
    top_protein = countAb(df_run_nodecoy_fixed.copy())
    top_protein_hybrid = {key:val for key,val in top_protein.items() if len(key.split("|")) != 3}
    native_list = list(set(top_protein.keys())-set(top_protein_hybrid.keys()))
    top_protein_native = {key:val for key,val in top_protein.items() if key in native_list and "Hybrid" not in key}
    ab_hybrid_df = pd.DataFrame(createAbDf(df_run_nodecoy_fixed.copy(), top_protein_hybrid, top_protein_native, df_hybrid.copy()))
    
    ab_hybrid_df['left+right'] = ab_hybrid_df.left_abundance.values+ab_hybrid_df.right_abundance.values
    ab_hybrid_df['left-right'] = np.absolute(ab_hybrid_df.left_abundance.values-ab_hybrid_df.right_abundance.values)
    ab_hybrid_df['sort'] = (ab_hybrid_df['left+right'].values-ab_hybrid_df['left-right'].values)#*ab_hybrid_df.hybrid.values
    ab_hybrid_df.sort_values(by=["sort","max_ions_matched","hybrid_abundance"], ascending=False, inplace=True)
    ab_hybrid_df = ab_hybrid_df.reset_index(drop=True)
    ab_hybrid_df.drop(columns=["sort","left+right","left-right"], inplace=True)
    ab_hybrid_df.to_csv(os.path.join(ut.checkDir(args.out),"hybrid_abundance_all.csv"))
    ab_hybrid_df = ab_hybrid_df[~ab_hybrid_df.name.str.contains("/")].reset_index(drop=True)
    ab_hybrid_df.to_csv(os.path.join(ut.checkDir(args.out),"hybrid_abundance.csv"))


if __name__ == '__main__':
    main()
