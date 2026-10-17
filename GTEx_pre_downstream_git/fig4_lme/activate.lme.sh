sbatch -o SMTSDbrain_lme.y.out -e SMTSDbrain_lme.y.err --wrap "Rscript lme.across_SMTSDbrain.R --inferred-covs yes"   # include InferredCov1/2
sbatch -o SMTSDbrain_lme.n.out -e SMTSDbrain_lme.n.err --wrap "Rscript lme.across_SMTSDbrain.R --inferred-covs no"    # exclude InferredCov1/2
sbatch -o lme.across_SMTS.out -e lme.across_SMTS.err --wrap "Rscript  lme.across_SMTS.R"
