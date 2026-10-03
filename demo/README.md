# Demo (demo folder) 
A standalone, self-contained notebook (independent of the SLURM pipeline scripts) that walks through the DIV diversity metric on both simulated and real GTEx-derived data:
- Simulated data — builds a small synthetic isoform-expression table and computes DIV_1 on it as a sanity check/illustration of the metric.
- Biological data — uses small demo subsets (first 10 samples) cut from the pipeline's real outputs (outputdir.1/readtpm_t.csv.gz for isoform-level TPM and outputdir.5/other_and_HK-RP/tpm.csv.gz for transcriptome-level TPM), with the exact zcat/cut commands used to generate them documented in the notebook.
1. **Isoform diversity** — computes per-gene isoform diversity (I_DIV) and visualizes, for an example sample, the isoform expression levels and relative abundances underlying the diversity score.
2. **Transcriptome diversity** — computes per-sample transcriptome-wide diversity (T_DIV) and plots, per sample, gene expression percentile vs. min-max-scaled transcriptome diversity.
It's meant as a lightweight, visual introduction to how the DIV metric (also used throughout the numbered pipeline scripts) behaves, rather than as part of the production pipeline itself.

**Expected output:** isoform and transcriptome diversity scores, visually described based on inform and gene expression respectively.
