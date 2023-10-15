#!/bin/bash
set -e # stop on error
set -u # raise error if variable is unset
set -o pipefail # fail if any prior step failed

echo "Program is running ..."
python ../../src/fix_protein.py --comet_in ../outputs/comet_run_2 --hypedsearch_in ../outputs/Hypedsearch_outputs --out ../outputs/comet_run_2
echo "The hybrid records in Comet output files have been fixed"
python ../../src/count_abundance.py --comet_in ../outputs/comet_run_2/fixed_comethypedsearch_outputs.csv --hypedsearch_in ../outputs/Hypedsearch_outputs --out ../outputs/results
echo "The hybrid abundance table has been saved"
python ../../src/plot_hybridAb.py --hybrid_in ../outputs/results/hybrid_abundance.csv --score max_ions_matched --outfile ../outputs/results/hybrid_Ab_ions.png
echo "The hybrid abundance figure has been saved"
