import pandas as pd
import numpy as np
import argparse
import warnings
import os
warnings.filterwarnings('ignore')


def quantile_normalize(df):
    # 1. Sort each column and record the rank of each value
    sorted_df = pd.DataFrame(np.sort(df.values, axis=0), columns=df.columns)
    
    # 2. Compute the row-wise mean of the sorted values (the "reference distribution")
    row_means = sorted_df.mean(axis=1)
    
    # 3. Replace each value with the mean corresponding to its rank
    rank_df = df.rank(method="min").astype(int) - 1  # 0-indexed ranks
    normalized = rank_df.apply(lambda col: col.map(row_means))
    
    return normalized

parser = argparse.ArgumentParser(description='Quantile normalization of dataframes')
parser.add_argument('--df', type=str, help='Dataframe')
args = parser.parse_args()


inputpath = os.getcwd() + '/outputdir.4/'
outputpath = os.getcwd() + '/outputdir.6/'

if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )

quantile_normalize(pd.read_csv(inputpath + args.df, index_col=0)).to_csv(outputpath + args.df, compression='gzip', index = True)