import subprocess
import sys
import os

# Read in the cell type pair interaction dataframe from each cohort, and pass each one as an argument to
# pheatmap.R. This Rscript will generate a heatmap for each assigned dataframe, and write it as a pdf file.

def main():
    cancer_type = sys.argv[1]
    pair_df_dir = "square_sum_results/results_" + cancer_type + "/interaction_tables"
    out_dir = "square_sum_results/results_" + cancer_type + "/interaction_figures"

    for file in os.listdir(pair_df_dir):
        print(file)
        cohort = file[9:len(file)].replace('.csv', '')
        print(cohort)
        r_cmd = ' '.join(['Rscript', 'generate_hmap.R', '-d', pair_df_dir, '-f', file, '-n', str(file[7]), '-c', cohort, '-o', out_dir])
        process = subprocess.Popen(r_cmd, shell=True).wait()

if __name__=="__main__":
    main()
