# input_files: public input data

Step 1 of the pipeline. All later steps read their raw input from this folder.

## Run

```bash
cd input_files
bash commands.sh
```

`commands.sh` downloads, into this folder:

| File | Source | Used by |
|------|--------|---------|
| `GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct.gz` | GTEx v8 portal | `GTEx_upstream` |
| `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` | GTEx v8 portal | `GTEx_upstream` |
| `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt` | GTEx v8 portal | `GTEx_upstream` |
| `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | GTEx v8 portal | `GTEx_upstream` |
| `gencode.v26.annotation.gff3` | GENCODE release 26 (the annotation of GTEx v8) | `GTEx_upstream` |
| `*.genes.results.gz`, `*.isoforms.results.gz` (from `GSE245325_RAW.tar`) | GEO GSE245325, Ulicevic et al. 2024 | `Ulicevic_upstream` |
| `gencode.v28.annotation.gff3` | GENCODE release 28 | `Ulicevic_upstream` |
| `isoform_txt.zip` | PRIDE PXD046084, Ulicevic et al. 2024 (proteomics) | `downstream/Ulicevic.ipynb` |

GTEx v8 open-access files require no login. `SubjectPhenotypesDS.txt` from the open-access portal holds
age as 10-year brackets, which is what the pipeline uses.

## Files included in this repository

- `ribosomal_genes.csv`: the 152 ribosome-related genes (Ensembl ID, gene name) used to split genes
  into ribosomal and non-ribosomal (T-DIV is computed on the non-ribosomal genes).
  They are the genes of the Enrichr library GO_Cellular_Component_2025 (human) in terms whose name
  contains "ribosom", mapped to Ensembl IDs with Ensembl Biomart. The list is saved because both
  resources change over time; `make_ribosomal_genes.py` shows how it was made.
- `Donovan2020/41467_2020_14561_MOESM{5..18}_ESM.csv`: cell-type composition estimates of GTEx samples,
  Supplementary Data 4–17 of Donovan et al. 2020, *Nat Commun* 11:955
  (https://www.nature.com/articles/s41467-020-14561-0), shared under CC BY 4.0. Used by
  `GTEx_upstream/5.calculate_transcriptome_div.py` and `downstream/transcriptome_div.ipynb`.

## File to add by hand

- `mele2015.AgeDiffExpGenes.S12.csv`: Table S12 of Melé et al. 2015, *Science* 348:660
  (age-associated genes), saved as CSV. Used only by `downstream/lme.ipynb` for the comparison with
  published age effects. *[ADD: how the CSV was made from the supplementary file: sheet, columns kept.]*
