import scanpy as sc
import os
import pandas as pd
import sys
import scipy
import numpy as np
from calculate_adj_bk3 import build_adj_df2
from calculate_adj_bk import build_adj_df
from cluster_analysis import find_clusters
from cell_type_interactions import interactions_by_clust, ctype_ratios
# from cluster_co_exprsn import find_clust_coexprsd



def main():
    # read list of h5ad-formatted single cell datasets, and convert to AnnData objects
    cancer_type = sys.argv[1]
    data_dir = "../sc_data/" + cancer_type 
    error_file = open("error_cohorts.txt", 'a')
   
    # file_lst = [os.listdir(data_dir)[1], os.listdir(data_dir)[3], os.listdir(data_dir)[0]]
    for file in os.listdir(data_dir):
        print(file)
        try:
            adata = sc.read_h5ad(os.path.join(data_dir, file))
        except ValueError:
            print(file + " has an error ")
            error_file.write(file + '\n')
            continue
        else:
            print(adata)
            cohort = file.replace('.h5ad', '')
            adj_csr, cell_names = build_adj_df2(adata)
            adj_df = pd.DataFrame(adj_csr.A, index= cell_names, columns=cell_names)
            print("after loading csr, the adj mat has shape " + str(adj_df.shape))
            scipy.sparse.save_npz("square_sum_results/results_" + cancer_type + "/adj_matrices/" + "adj_csr_" + cohort + ".npz", adj_csr)
            clust_df, cell_annot, exprsn_filt = find_clusters(adata, adj_df)
            clust_df.to_csv("square_sum_results/results_" + cancer_type + "/cluster_assignments/" + "clust_df_" + cohort + ".csv")
            cell_annot.to_csv("square_sum_results/results_" + cancer_type + "/cell_annotations/" + "cell_annot_" + cohort + ".csv")
            exprsn_filt.to_csv("square_sum_results/results_" + cancer_type + "/expression_matrices/" + "exprsn_filt_" + cohort + ".csv")
            print("starting cell pair analysis")
            pair_df1, pair_df2 = interactions_by_clust(adj_df, cell_annot, clust_df)
            pair_df1.to_csv("square_sum_results/results_" + cancer_type + "/interaction_tables/pair_df1_" + cohort + ".csv")
            pair_df2.to_csv("square_sum_results/results_" + cancer_type + "/interaction_tables/pair_df2_" + cohort + ".csv")
            print("starting cluster composition")
            plt = ctype_ratios(cell_annot, clust_df)
            plt.savefig("square_sum_results/results_" + cancer_type + "/cell_type_ratios/clust_composition_plt_" + cohort + ".pdf")
        
        
    error_file.close()

if __name__=="__main__":
    main()