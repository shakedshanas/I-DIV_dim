# demo: isoform and transcriptome diversity on a small example

A self-contained example of the diversity index of the paper. It needs only the Python environment of
`../environment.yml` (numpy, pandas, matplotlib, seaborn, scikit-learn, jupyter) and the two CSV files in
this folder; `DIV()` is defined inside the notebook, identical to `../div.py`.

## Data

The first 9 GTEx v8 samples (all from donor GTEX-111FC):

- `input_for_isoform_div.csv`: transcript TPM; one row per transcript, `accession_number` = its gene.
- `input_for_transcriptome_div.csv`: gene TPM; one row per gene.

## Run

```bash
cd demo
jupyter nbconvert --to notebook --execute demo.ipynb --output demo.executed.ipynb
```

or open `demo.ipynb` in Jupyter and *Run All*. Run time: 41 seconds on a login node of our cluster.

## Expected output

1. **Simulated data:** 50 vectors of 10 isoforms whose expression becomes more and more dominated by one
   isoform; their I-DIV falls from 0.78 to 0.61, plotted as a line.
2. **I-DIV:** `I_DIV`, a genes × samples table of isoform diversity, and bar plots of the isoform abundances
   of three example genes (3, 5 and 10 expressed isoforms) in sample GTEX-111FC-0126-SM-5N9DL.
3. **T-DIV:** one value per sample:

   ```
   GTEX-111FC-0126-SM-5N9DL    0.743522
   GTEX-111FC-0226-SM-5N9B8    0.864837
   GTEX-111FC-0326-SM-5GZZ1    0.634731
   GTEX-111FC-0426-SM-5N9CV    0.809148
   GTEX-111FC-0526-SM-5GZZ8    0.870944
   GTEX-111FC-0626-SM-5N9CU    0.788239
   GTEX-111FC-0826-SM-5GZWO    0.783891
   GTEX-111FC-1126-SM-5GZWU    0.840615
   GTEX-111FC-2526-SM-5GZXU    0.862663
   ```

   and, per sample, the normalised gene expression plotted against its percentile.

The committed `demo.ipynb` contains these outputs; a fresh run reproduces the T-DIV values exactly.
