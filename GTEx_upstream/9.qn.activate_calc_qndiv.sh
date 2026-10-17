#!/bin/bash
#SBATCH --job-name=div
#SBATCH --array=0-10000
#SBATCH --partition=hive7d
#SBATCH --output=outputs/%A_%a.DIV.output
# run after 7.qn.activate_jointables.sh and 8.make_headers.sh have finished

ngenes=$(zcat 4.outputs/gene_annotation.csv.gz | wc -l)
smallidx=$(($ngenes - 10000))

# for df in QN.Brain.I.csv.gz QN.SMTS.I.csv.gz QN.SMTSD.G.csv.gz; do
#     python 9.qn.calc_qndiv.py --geneidx $SLURM_ARRAY_TASK_ID --iinput $df
#     if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
#         ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))
#         python 9.qn.calc_qndiv.py --geneidx $ACTUAL_JOB_ID --iinput $df

#     fi
# done


python 9.qn.calc_qndiv.py --geneidx $SLURM_ARRAY_TASK_ID --iinput QN.SMTS.I.csv.gz
python 9.qn.calc_qndiv.py --geneidx $SLURM_ARRAY_TASK_ID --iinput QN.SMTSD.I.csv.gz
python 9.qn.calc_qndiv.py --geneidx $SLURM_ARRAY_TASK_ID --iinput QN.Brain.I.csv.gz
if [ $SLURM_ARRAY_TASK_ID -le $smallidx ] ; then 
    ACTUAL_JOB_ID=$((SLURM_ARRAY_TASK_ID + 10000))
    python 9.qn.calc_qndiv.py --geneidx $ACTUAL_JOB_ID --iinput QN.SMTS.I.csv.gz
    python 9.qn.calc_qndiv.py --geneidx $ACTUAL_JOB_ID --iinput QN.SMTSD.I.csv.gz
    python 9.qn.calc_qndiv.py --geneidx $ACTUAL_JOB_ID --iinput QN.Brain.I.csv.gz

fi
