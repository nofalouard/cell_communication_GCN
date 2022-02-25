import scanpy as sc
import os
import pandas as pd
# from calculate_adj_bk import build_adj_df
# from cluster_analysis import find_clusters
# from cell_type_interactions import identify_interactions
from positive_control_ranks import lig_rec_ranks
from matplotlib import pyplot as plt



def main():
    # read list of h5ad-formatted single cell datasets, and convert to AnnData objects
    data_dir = "../sc_data" 
    adata_lst = []
    file_lst = [os.listdir(data_dir)[2], os.listdir(data_dir)[1], os.listdir(data_dir)[0]]
    # file_lst = ['NSCLC_GSE153935.h5ad', 'NSCLC_GSE150660.h5ad', 'NSCLC_GSE143423.h5ad', 'NSCLC_GSE117570.h5ad', 'NSCLC_GSE99254.h5ad']
    for file in file_lst[1:3]:
        print(file)
        adata = sc.read_h5ad(os.path.join(data_dir, file))
        print(adata)
        cohort = file.replace('.h5ad', '')
        # adj_df = pd.read_csv("results/adj_df_" + cohort + ".csv")
        # print("finished reading the adj_df")
        # cell_annot = pd.read_csv("results/cell_annot_" + cohort + ".csv")
        # exprsn_filt = pd.read_csv("exprsn_filt_" + cohort + ".csv")
        # pair_df = identify_interactions(adj_df, cell_annot)
        # pair_df.to_csv("pair_df_" + cohort + ".csv")
        plt = lig_rec_ranks(adata)
        plt.savefig("pos_control_ranks_" + cohort + ".pdf")


if __name__=="__main__":
    main()