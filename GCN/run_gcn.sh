#!/bin/bash
#SBATCH --job-name=GCN
#SBATCH --mem=96G       # total memory need
#SBATCH --time=96:00:00
#SBATCH -c 6 #number of cores
#SBATCH -o o_%j.txt
#SBATCH -e e_%j.txt


source activate gcn
all_dirs=(../sc_data/*)
for dir in "${all_dirs[@]:31}"
do 
    dir=${dir%*/} 
    dir="${dir##*/}" 
    echo $dir
    python analyze_adj.py $dir
    echo "finished adj calculations"
    python cell_pair_hmap.py $dir
    echo "finished cell pair analysis"
    python analyze_coexpression.py $dir
    echo "finished coexpression analysis"
done