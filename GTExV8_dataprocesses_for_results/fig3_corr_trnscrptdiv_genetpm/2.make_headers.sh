#!/bin/bash

python T-DIV.py -i 0 -s all   -f other_and_HK-RP -z rows
python T-DIV.py -i 0 -s SMTSD   -f other_and_HK-RP -z rows
python T-DIV.py -i 0 -s SMTS   -f other_and_HK-RP -z rows
python mI-DIV.py -i 0 -s all -z rows
python mI-DIV.py -i 0 -s SMTSD -z rows
python mI-DIV.py -i 0 -s SMTS  -z rows