#!/bin/bash

# sbatch -J SMTSD -o outputs/joindfs.SMTSD.IDIV.out -e outputs/joindfs.SMTSD.IDIV.err --wrap "python joindfs.py --data IDIV --basedon SMTSD"
sbatch -J SMTSD -o outputs/joindfs.SMTSD.G.out -e outputs/joindfs.SMTSD.G.err --wrap "python  7.qn.joindfs.py --data G --basedon SMTSD"
sbatch -J SMTSD -o outputs/joindfs.SMTSD.I.out -e outputs/joindfs.SMTSD.I.err --wrap "python  7.qn.joindfs.py --data I --basedon SMTSD"

sbatch -J SMTSDbrain -o outputs/joindfs.SMTSDbrain.G.out -e outputs/joindfs.SMTSDbrain.G.err --wrap "python  7.qn.joindfs.py --data G --basedon SMTSDbrain"
sbatch -J SMTSDbrain -o outputs/joindfs.SMTSDbrain.I.out -e outputs/joindfs.SMTSDbrain.I.err --wrap "python  7.qn.joindfs.py --data I --basedon SMTSDbrain"

# sbatch -J SMTS -o outputs/joindfs.SMTS.IDIV.out -e outputs/joindfs.SMTS.IDIV.err --wrap "python joindfs.py --data IDIV --basedon SMTS"
sbatch -J SMTS -o outputs/joindfs.SMTS.G.out -e outputs/joindfs.SMTS.G.err --wrap "python  7.qn.joindfs.py --data G --basedon SMTS"
sbatch -J SMTS -o outputs/joindfs.SMTS.I.out -e outputs/joindfs.SMTS.I.err --wrap "python  7.qn.joindfs.py  --data I --basedon SMTS"