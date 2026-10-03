import pandas as pd
import numpy as np
import scipy.stats as st



tpm = pd.read_csv(f'/lustre1/home/martinmikl/sshanas/outsourcedata/pipeline_2025/outputdir.3/readtpm_g.csv.gz', compression='gzip', index_col=0)
tpm = np.log2(tpm + 1)
# z-score by sample
tpm.apply(st.zscore, axis = 'rows').to_csv('readtpm_g.zscore.rows.csv.gz', index=True, compression='gzip') # every sample: mean 0 st 1
# tpm.apply(st.zscore, axis = 'columns').to_csv('readtpm_g.zscore.columns.csv.gz', index=True, compression='gzip') # every gene: mean 0 st 1