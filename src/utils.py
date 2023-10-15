import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt

def checkDir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
    return dir_path


def getData(dir_path, skiprow1=True):
    all_files = glob(os.path.join(dir_path, "*.txt"))
    if skiprow1 == True:
        df = pd.concat((pd.read_csv(f,sep='\t',skiprows=1) for f in all_files), ignore_index=True)
    else:
        df = pd.concat((pd.read_csv(f,sep='\t') for f in all_files), ignore_index=True)
    return df


def dropDup(x):
  return list(dict.fromkeys(x))


def checkSeq(seq, max_len, pre_aa, next_aa, output_df):
  result = ""
  if len(seq.split("-")[0])+len(seq.split("-")[1]) <= max_len:
    if seq.split("-")[0][0] == pre_aa:
      if next_aa in output_df["next aa"][output_df.sequence == seq].values[0]:
        result = seq
  return result


def getNativeAb(side, dict_ab):
  ab_list = [val for key,val in dict_ab.items() if key.split('|')[-1] in side]
  if len(ab_list) != 0:
    return max(ab_list)
  else:
    return 0
