import scanpy as sc
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def lig_rec_ranks(adata):
    exprsn = pd.DataFrame(adata.X.A, index = adata.obs.index, columns = adata.var_names)
    
    # read in list of proteins that serve as receptors and ligands, and file identifying which pairs interact
    prot_file = pd.read_csv("../interaction_data/protein_input.csv")
    gene_file = pd.read_csv("../interaction_data/gene_input.csv")
    inter_file = pd.read_csv("../interaction_data/interaction_input.csv")
    lig_list = list(prot_file[prot_file["receptor"]==False]["uniprot"])
    rec_list = list(prot_file[prot_file["receptor"]==True]["uniprot"])

    simple_inter = inter_file[inter_file['partner_a'].isin(list(prot_file["uniprot"])) & \
                   inter_file['partner_b'].isin(list(prot_file["uniprot"]))]
    parta = list(simple_inter['partner_a'])
    partb = list(simple_inter['partner_b'])
    parta_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), parta))

    lig_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), lig_list))
    rec_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), rec_list))

    lig_df = pd.DataFrame(list(zip(lig_list, lig_genes)))
    rec_df = pd.DataFrame(list(zip(rec_list, rec_genes)))

    gene_df = pd.DataFrame(list(zip(parta, parta_genes)))

    partb_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), partb))
    gene_df2 = pd.DataFrame(list(zip(partb, partb_genes)))

    gene_df[1]=gene_df[1].apply(lambda row: gene_file.loc[row, 'gene_name'].to_list()[0])

    gene_df2[1]=gene_df2[1].apply(lambda row: gene_file.loc[row, 'gene_name'].to_list()[0])

    simple_inter.loc[:,'parta_gene'] = gene_df.loc[:,1]

    simple_inter.loc[:,'partb_gene'] = gene_df2.loc[:,1]
    parta = [x for x in gene_df[1] if x in exprsn.columns]
    partb = [x for x in gene_df2[1] if x in exprsn.columns]

    exprsn_parta = exprsn[parta]
    exprsn_partb = exprsn[partb]

    exprsn_filt = pd.concat([exprsn_parta,exprsn_partb], axis=1)

    exprsn_filt=exprsn_filt.T.drop_duplicates().T

    all_pairs = []
    partners = zip(list(gene_df[1]), list(gene_df2[1]))
    for a,b in partners:
        all_pairs.append((a,b))

    pairs_found = []
    for pair in all_pairs:
        if pair[0] in exprsn_filt.columns and pair[1] in exprsn_filt.columns:
            pairs_found.append(pair)

    gene_lst = list(exprsn.columns)

    num_genes = exprsn.shape[1]

    # Identify each gene pair's corresponding index in the flattened gene correlation matrix. 
    # Method: For each pair, find partner a's index in the gene list, and multiply by the total number of genes.
    #         Then, add this result to partner b's index in the gene list.

    flat_ind = []
    for pair in pairs_found:
        inda = gene_lst.index(pair[0])
        indb = gene_lst.index(pair[1])
        flat_ind.append(inda*num_genes + indb)
    
    exprsn_mat = adata.X.A
    corr_mat = np.corrcoef(exprsn_mat, rowvar=False)
    sorted_corr = np.argsort(corr_mat, axis = None).tolist()
    # For each ligand-receptor pair found in the expression matrix, identify where its correlation coefficient is ranked
    # in comparison to the sorted coefficients (increasing order) of the other genes. Save rank values in a dictionary.
    print(len(flat_ind))
    pair_ranks = {}
    total_corrs = len(sorted_corr)
    for num, ind in enumerate(flat_ind):
        rank = sorted_corr.index(ind)
        percent = rank/total_corrs*100
        pair_ranks[pairs_found[num]] = percent
    print('making plot')
    ranks_df = pd.Series(list(pair_ranks.values()), index = list(pair_ranks.keys()))
    n_bins = 20
    fig, axs = plt.subplots(1, 1,figsize =(10, 7),tight_layout = True)
    plt.xlabel("Percentile")
    plt.ylabel("Frequency")
    plt.title('Percentiles of Ligand-receptor Correlation Coefficients')
    plt.xticks(np.arange(0, 110, step=10))
    axs.hist(ranks_df.values, bins = n_bins)
    return plt
    



    



    












