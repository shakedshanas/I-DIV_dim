import pandas as pd
import numpy as np
import argparse
from sklearn.linear_model import LinearRegression
import os
import warnings
warnings.filterwarnings('ignore')

inputpath = os.getcwd() + '/outputdir.3/'
outputpath = os.getcwd() + '/outputdir.4/'

if not os.path.exists( outputpath ):
    os.mkdir(outputpath  )

parser = argparse.ArgumentParser(description='Create diversity files')
parser.add_argument('--dataname', type=str, choices=['readtpm_g', 'div1'], help="name of the data to be corrected,")
parser.add_argument('--grouptype', type=str, choices=['i', 'ii'], help="i corrects the gene for all the samples ; ii corrects per sample group")
parser.add_argument('--basedon', type=str, choices=['SMTS', 'SMTSD'], help="if used ii, will it be basedon SMTS or SMTSD")
parser.add_argument('--lgtype', type=str, choices=['1', '2'], help="2 is correcting only by the technical properties ; 1 is correcting and putting back")
# 'SMTS': 'tissue', 'SMTSD': 'subtissue', 'SMNABTCHT': 'isolation', 'SMTSISCH': 'ischemic', 'SMRIN': 'smpl_quality', 'SMGEBTCH': 'batch_id', 'SEX': 'sex', 'AGE': 'age', 'DTHHRDY': 'death'
args = parser.parse_args()

metadata = pd.read_csv('outputdir.1/sample_metadata.csv.gz', compression='gzip', index_col=['SAMPID'])
metadata = metadata.assign(deathtype2 = metadata.DTHHRDY.apply(lambda x: x ==2).astype(int))
metadata = metadata.assign(female = metadata.SEX.apply(lambda x: x ==2).astype(int))

tissue_matrix_preReg = pd.get_dummies(metadata['SMGEBTCH']).astype(int)

if args.lgtype == str(1):
    tissue_matrix_preReg = tissue_matrix_preReg.join(metadata[['SMRIN', 'SMTSISCH', 'age_decade', 'deathtype2', 'female']]) 
    # tissue_matrix_preReg = tissue_matrix_preReg.join(metadata[['SMRIN', 'SMTSISCH', 'deathtype2', 'female']])
elif args.lgtype == str(2):
    tissue_matrix_preReg = tissue_matrix_preReg.join(metadata[['SMRIN', 'SMTSISCH', 'deathtype2']])

def load_matrix(samplelist):
    # p = '/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/'
    if args.dataname == 'div1':
        gtex = pd.read_csv(inputpath + 'div1.csv.gz', index_col=['Unnamed: 0'], usecols = ['Unnamed: 0'] + samplelist)
    if args.dataname == 'readtpm_g':
        gtex = pd.read_csv(inputpath + 'readtpm_g.csv.gz', index_col=['accession_number'], usecols = ['accession_number'] + samplelist)
        gtex = np.log2(gtex + 1)
    return gtex


reg = LinearRegression()
# %%
    
def lreg_2(samplelist=None):
    
    resid_mat = pd.DataFrame(index=gtex.columns)
    predd_mat = pd.DataFrame(index=gtex.columns)
    temp_tissue_matrix_preReg = tissue_matrix_preReg.loc[samplelist]
    X = temp_tissue_matrix_preReg.values
    Y = gtex.values.T
    reg.fit(X,Y)

    # Get predictions and residuals
    predictions = reg.predict(X)  # (n_samples, n_genes)
    residuals = Y - predictions

    # Convert to DataFrames with proper indexing
    resid_mat = pd.DataFrame(residuals.T, index=gtex.index, columns=gtex.columns)
    predd_mat = pd.DataFrame(predictions.T, index=gtex.index, columns=gtex.columns)
    
    return resid_mat, predd_mat



def lreg_1(correct_vars=['age_decade', 'female'], samplelist=None):
# def lreg_1(correct_vars=['female'], samplelist=None):
    """
    Optimized loop version - much faster than original while keeping loop structure
    """
    temp_tissue_matrix_preReg = tissue_matrix_preReg.loc[samplelist]
    resid_mat = pd.DataFrame(index=gtex.columns)
    predd_mat = pd.DataFrame(index=gtex.columns)
    
    # Pre-allocate output matrices
    n_genes, n_samples = gtex.shape
    n_features = temp_tissue_matrix_preReg.shape[1]
    
    residuals = np.zeros((n_genes, n_samples))
    predictions = np.zeros((n_genes, n_samples))
    coefficients = np.zeros((n_genes, n_features))
    
    # Convert to numpy once (avoid repeated conversions)
    X = temp_tissue_matrix_preReg.values
    feature_names = temp_tissue_matrix_preReg.columns.tolist()
    
    # Find indices of correction variables
    correct_indices = [feature_names.index(var) for var in correct_vars]
    
    # Single regression object
    reg = LinearRegression()
    
    # Loop through genes
    for i, gene in enumerate(gtex.index):
        y = gtex.loc[gene].values
        
        # Fit regression
        reg.fit(X, y)
        pred = reg.predict(X)
        resid = y - pred
        
        # Add back correction variables
        for var_idx in correct_indices:
            resid += X[:, var_idx] * reg.coef_[var_idx]
        
        # Store results
        residuals[i] = resid
        predictions[i] = pred
        coefficients[i] = reg.coef_
    
    # Convert to DataFrames
    resid_mat = pd.DataFrame(residuals, index=gtex.index, columns=gtex.columns)
    predd_mat = pd.DataFrame(predictions, index=gtex.index, columns=gtex.columns)
    # coef_df = pd.DataFrame(coefficients, index=gtex.index, columns=feature_names)
    
    return resid_mat, predd_mat # , coef_df

if args.grouptype == 'i':
    l = metadata.index.tolist()
    gtex = load_matrix(samplelist=l)
    resid_mat, predd_mat = globals()['lreg_' + str(args.lgtype)](samplelist=l)
elif args.grouptype == 'ii':
    for subt in metadata[args.basedon].unique():
        l = metadata[metadata[args.basedon] == subt].index.tolist()
        gtex = load_matrix(samplelist=l)
        resid_mat, predd_mat = globals()['lreg_' + str(args.lgtype)](samplelist=l)
        
        for ffeature in ['resid', 'predd']:
            filename = '.'.join([args.grouptype, str(args.lgtype), args.basedon, ffeature, args.dataname.split('.')[0], 'csv.gz'])
            if os.path.exists(outputpath + filename):
                # Read existing file and join with new data
                existing_df = pd.read_csv(outputpath + filename,index_col=0, compression='gzip')
                combined_df = existing_df.join(globals()[ffeature + '_mat'])
                combined_df.to_csv(outputpath + filename, index=True, compression='gzip')
            else:
                # Create new file
                globals()[ffeature + '_mat'].to_csv(outputpath + filename, index=True, compression='gzip')

# %%
if args.grouptype == 'i':
    resid_mat.to_csv(outputpath + '.'.join([args.grouptype, args.lgtype, 'resid', args.dataname.split('.')[0], 'csv.gz']), index=True, compression='gzip')
    predd_mat.to_csv(outputpath + '.'.join([args.grouptype, args.lgtype, 'pred', args.dataname.split('.')[0], 'csv.gz']), index=True, compression='gzip')
# elif args.grouptype == 'ii':
#     resid_mat.to_csv('.'.join([args.grouptype, args.lgtype, args.basedon, 'resid', args.dataname.split('.')[0], 'csv.gz']), index=True, compression='gzip')
#     predd_mat.to_csv('.'.join([args.grouptype, args.lgtype, args.basedon, 'pred', args.dataname.split('.')[0], 'csv.gz']), index=True, compression='gzip')
        
