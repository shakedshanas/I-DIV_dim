#!/bin/bash
# three linear mixed-effects models (random intercept per donor), one job each
sbatch -o lme.across_SMTSDbrain.out -e lme.across_SMTSDbrain.err --wrap "Rscript lme.across_SMTSDbrain.R"
sbatch -o lme.across_SMTS.out -e lme.across_SMTS.err --wrap "Rscript lme.across_SMTS.R"
sbatch -o lme_across_SMTS_exRINISCH.out -e lme_across_SMTS_exRINISCH.err --wrap "Rscript lme_across_SMTS_exRINISCH.R"
