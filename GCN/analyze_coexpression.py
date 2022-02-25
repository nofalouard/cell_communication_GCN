import scanpy as sc
import os
import pandas as pd
import sys
from calculate_adj_bk import build_adj_df
# from cluster_analysis import find_clusters
# from cell_type_interactions import identify_interactions
from cluster_co_exprsn import find_clust_coexprsd, find_pair_coexprsd



def main():
    cancer_type = sys.argv[1]
    data_dir = "../sc_data/" + cancer_type
    # read list of h5ad-formatted single cell datasets, and convert to AnnData objects
    
    clust_coexp = []
    for file in os.listdir(data_dir):
        print(file)
        try:
            adata = sc.read_h5ad(os.path.join(data_dir, file))
        except ValueError:
            print(file + " has an error ")
            continue
        else:
            cohort = file.replace('.h5ad', '')
            clust_df = pd.read_csv("square_sum_results/results_" + cancer_type + "/cluster_assignments/" + "clust_df_" + cohort + ".csv")
            cohort_coexprsd = find_clust_coexprsd(clust_df, adata)
            # pair = tuple(['M1', 'CD4Tnaive'])
            # pair_coexprsd, pair_found = find_pair_coexprsd(clust_df, clust_annot, pair, adata)
            # if pair_found == True:
                # if melanoma_cpair_coexp != []:
                    # melanoma_cpair_coexp = list(set(pair_coexprsd) & set(melanoma_cpair_coexp))
                # else:
                    # melanoma_cpair_coexp = pair_coexprsd 
            # else:
                # print("cohort doesnt have desired cell type pair")
            print("there are " + str(len(cohort_coexprsd)) + " cohort co-expressed gene pairs ")
            if clust_coexp != []:
                clust_coexp = list(set(cohort_coexprsd) & set(clust_coexp))
            else:
                clust_coexp = cohort_coexprsd 
            print( str(len(clust_coexp)) + " are shared among all cohorts so far")
    print("across all " + cancer_type +  " cohorts, there are " + str(len(clust_coexp)) + " shared gene pairs")
    pair_df = pd.DataFrame(clust_coexp)
    pair_df.to_csv("square_sum_results/results_" + cancer_type + "/co_expression/" + "coexprsd_genes_" + cancer_type + ".csv")


if __name__=="__main__":
    main()