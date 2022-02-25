import scanpy as sc
#from anndata import AnnData
import pandas as pd
import numpy as np
import scipy.stats as stats
import scipy
from itertools import compress


# retrieve normalized gene exprsn matrix
# convert matrix to dataframe with columns as gene names and row indices as cells  
# adata = sc.read_h5ad("../data/skcm_mat.h5ad")
# sc.pp.normalize_total(adata)
def build_adj_df2(adata):
    
    exprsn = pd.DataFrame(adata.X.A)
    exprsn.set_index(adata.obs.index, inplace=True)
    exprsn.columns = adata.var_names
    
    # from protein and gene lists, mark which gene pairs code for interacting ligand-protein
    # complexes
    prot_file = pd.read_csv("../interaction_data/protein_input.csv")
    gene_file = pd.read_csv("../interaction_data/gene_input.csv")
    inter_file = pd.read_csv("../interaction_data/interaction_input.csv")
    lig_list = list(prot_file[prot_file["receptor"]==False]["uniprot"])
    rec_list = list(prot_file[prot_file["receptor"]==True]["uniprot"])
    simple_inter = inter_file[inter_file['partner_a'].isin(list(prot_file["uniprot"])) & inter_file['partner_b'].isin(list(prot_file["uniprot"])) ]
    parta = list(simple_inter['partner_a'])
    partb = list(simple_inter['partner_b'])
    
    parta_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), parta))
    partb_genes = list(map(lambda gene: gene_file.index[gene_file['uniprot']==gene].to_list(), partb))
    gene_df = pd.DataFrame(list(zip(parta, parta_genes)))
    gene_df2 = pd.DataFrame(list(zip(partb, partb_genes)))
    gene_df[1]=gene_df[1].apply(lambda row: gene_file.loc[row, 'gene_name'].to_list()[0])
    gene_df2[1]=gene_df2[1].apply(lambda row: gene_file.loc[row, 'gene_name'].to_list()[0])

    
    simple_inter.loc[:,'parta_gene'] = gene_df.loc[:,1]
    simple_inter.loc[:,'partb_gene'] = gene_df2.loc[:,1]
    parta_lost = [x for x in gene_df[1] if x not in exprsn.columns]
    partb_lost = [x for x in gene_df2[1] if x not in exprsn.columns]
    parta = [x for x in gene_df[1] if x in exprsn.columns]
    partb = [x for x in gene_df2[1] if x in exprsn.columns]
    

    exprsn_parta = exprsn[parta]
    exprsn_partb = exprsn[partb]


    # in expression matrix, only keep genes that code for any of the marked ligand-receptor pairs 
    exprsn_filt = pd.concat([exprsn_parta,exprsn_partb], axis=1)
    exprsn_filt=exprsn_filt.T.drop_duplicates().T
    
    

    # generate 
    all_pairs = []
    partners = zip(list(gene_df[1]), list(gene_df2[1]))
    for a,b in partners:
        all_pairs.append((a,b))

    pairs_found = []
    for pair in all_pairs:
        if pair[0] in exprsn_filt.columns and pair[1] in exprsn_filt.columns:
            pairs_found.append(pair)


    # create 2 dataframes, one listing each cell's exprsn of the first element of each gene pair, and one for 
    # each cell's exprsn of the second element of each pair
    pair_exprsn1 = pd.DataFrame(columns=pairs_found, index = exprsn.index)
    pair_exprsn2 = pd.DataFrame(columns=pairs_found, index = exprsn.index)

    # # print(pair_exprsn1.shape, pair_exprsn2.shape)

    genes1 = [x[0] for x in pair_exprsn1.columns]
    for index, row in pair_exprsn1.iterrows():
        gene1_exprsn = list(map(lambda g: exprsn_filt.loc[index, g], genes1))
        pair_exprsn1.loc[index]= gene1_exprsn


    genes2 = [x[1] for x in pair_exprsn2.columns]
    for index, row in pair_exprsn2.iterrows():
        gene2_exprsn = list(map(lambda g: exprsn_filt.loc[index, g], genes2))
        pair_exprsn2.loc[index]= gene2_exprsn


    # construct dataframe with cell names as rows and columns. each element will store, in list format, a given cell pair's joint exprsn 
    # of all lig-rec pairs listed in the pair_exprsn dataframe columns
    # This dataframe will be converted to a 3d matrix, to be later used to build adjacency matrix

    num_cells = adata.X.A.shape[0]
    adj_mat = np.zeros((num_cells, num_cells))
    all_cells = list(exprsn.index) 
    print("here are num of cells")
    print(pair_exprsn2.shape[0])
    thresh = round(pair_exprsn2.shape[0]*0.05)
    # print(pair_exprsn2.shape[0])
    print(thresh)
    num_pairs = pair_exprsn1.shape[1]
    print("here are the num of pairs found")
    print(num_pairs)
    all_pair_scores = list(map(lambda x, y: gene_pair_score(x,adj_mat, pair_exprsn1, pair_exprsn2, all_cells), range(num_pairs),adj_mat ))

    # for lr_index in range(pair_exprsn2.shape[1]):
    #     p1_mean = pair_exprsn1.iloc[:,lr_index].mean()
    #     p2_mean = pair_exprsn2.iloc[:,lr_index].mean()
    #     dfa = pair_exprsn1[pair_exprsn1.iloc[:,lr_index] >p1_mean].iloc[:,lr_index]
    #     dfb = pair_exprsn2[pair_exprsn2.iloc[:,lr_index] >p2_mean].iloc[:,lr_index]
    #     # thresh = round(len(pair_exprsn1.iloc[:,lr_index])*0.1)
    #     # thresh = 2000
    #     # print(len(dfa), len(dfb))
    #     # print(thresh)
    #     if len(dfa) > thresh:
    #         dfa = dfa.astype(float).nlargest(thresh)
    #     if len(dfb) > thresh:
    #         dfb = dfb.astype(float).nlargest(thresh)
    #     if len(dfa) < len(dfb):
    #         adj_mat = pair_sum_over_dfa(adj_mat, dfa, dfb, all_cells)
    #     else:
    #         adj_mat = pair_sum_over_dfb(adj_mat, dfa, dfb, all_cells)
    print("finished with adj calculations")
    cell_inds = dict((value, idx) for idx,value in enumerate(all_cells))
    # adj_df = pd.DataFrame(adj_mat, index = all_cells, columns = all_cells)
    print("here")
    row_empty =list(~np.all(adj_mat == 0, axis=1))
    col_empty =list(~np.all(adj_mat == 0, axis=0))
    print("got empty rows/cols")
    # row_empty = list((~(adj_df == 0).all(axis=1)).values)
    # col_empty = list((~(adj_df == 0).all(axis=0)).values)

    rows_keep = list(compress(all_cells, row_empty))
    cols_keep = list(compress(all_cells, col_empty))
    cells_keep = list(set(rows_keep) & set(cols_keep))
    inds_keep = [cell_inds[x] for x in cells_keep]
    print("filtering adj mat")
    # adj_csr = scipy.sparse.csr_matrix(adj_mat)
    adj_mat = adj_mat[inds_keep,:][:, inds_keep]
    # adj_mat = adj_csr[inds_keep,:][:, inds_keep].A
    # adj_df = pd.DataFrame(adj_mat, index = cells_keep, columns = cells_keep)
    # print("here2")
    # # adj_df = adj_df[~empty_rows]
    # print("here3")
    # # adj_df = adj_df.loc[:,~empty_rows]
    # print("here4")
    # empty_cols = list((adj_df == 0).all(axis=0).index)
    # empty_cols_mat = [cell_inds[x] for x in empty_cols]
    # print("here5")
    # # adj_df = adj_df[~empty_cols]
    # print("here5")
    # adj_df = adj_df.loc[:,~empty_cols]
    print("after filtering, adj df has shape " + str(adj_mat.shape))
    # adj_csr = scipy.sparse.csr_matrix(adj_df.to_numpy())
    adj_csr = scipy.sparse.csr_matrix(adj_mat)
    return adj_csr, cells_keep
    # return adj_df

    
    
def geom_mean(x, y):
    return (x * y) ** 0.5

def square_sum(x, y):
    return (x**2 + y **2) ** 0.5

# Input: an index corresponding to one of the ligand receptor gene pairs
# For each cell pair combination, function calculates the distance formula result on the first cell's ligand expression
# and second cell's receptor expression
# Returns: a cell x cell pandas dataframe storing all cell pair distances associated with the inputted ligand-receptor pair 


def cell_cell_exprsn(df, lr_index, pair_exprsn1, pair_exprsn2):
    p1_mean = pair_exprsn1.iloc[:,lr_index].mean()
    p2_mean = pair_exprsn2.iloc[:,lr_index].mean()
    dfa = pair_exprsn1[pair_exprsn1.iloc[:,lr_index] >p1_mean].iloc[:,lr_index]
    dfb = pair_exprsn2[pair_exprsn2.iloc[:,lr_index] >p2_mean].iloc[:,lr_index]
    thresh = round(len(pair_exprsn1.iloc[:,lr_index])*0.01)
    
    # thresh = 2500
    if len(dfa) > thresh:
        dfa = dfa.astype(float).nlargest(thresh)
    if len(dfb) > thresh:
        dfb = dfb.astype(float).nlargest(thresh)
    df_sub = pd.DataFrame(columns=dfa.index, index = dfb.index)
    # print(df_sub.shape)
    df_sub=df_sub.apply(lambda row: pd.DataFrame(row).apply(lambda x: square_sum(dfa.loc[row.name], dfb.loc[x.name]),axis=1))
    # cell_cell = sum([df_sub, cell_cell])
    df_sub = df_sub.add(df, fill_value=0)
    return df_sub

def pair_sum_over_dfa(adj_mat, dfa, dfb, all_cells):
    
    inds_a = np.where(np.isin(all_cells,list(dfa.index)))[0]
    inds_b = np.where(np.isin(all_cells,list(dfb.index)))[0]
    shuff_inds = list(inds_b)
    dfb_vals = list(dfb.values)
    started = False
    while np.array_equal(shuff_inds, inds_b) == False or started == False:
        started = True
        sum_pair = square_sum(np.array(dfb_vals[0:len(dfa)]), dfa.values).astype(float)
        adj_mat[inds_a, shuff_inds[0:len(dfa)]] += sum_pair
        shuff_inds = shuff_inds[-1:] + shuff_inds[:-1]
        dfb_vals = dfb_vals[-1:] + dfb_vals[:-1]
    return adj_mat

def pair_sum_over_dfb(adj_mat, dfa, dfb, all_cells):

    inds_a = np.where(np.isin(all_cells,list(dfa.index)))[0]
    inds_b = np.where(np.isin(all_cells,list(dfb.index)))[0]
    shuff_inds = list(inds_a)
    dfa_vals = list(dfa.values)
    started = False
    while np.array_equal(shuff_inds, inds_a) == False or started == False:
        started = True
        sum_pair = square_sum(np.array(dfa_vals[0:len(dfb)]), dfb.values).astype(float)
        adj_mat[shuff_inds[0:len(dfb)], inds_b] += sum_pair
        shuff_inds = shuff_inds[-1:] + shuff_inds[:-1]
        dfa_vals = dfa_vals[-1:] + dfa_vals[:-1]
    return adj_mat

def gene_pair_score(lr_index, adj_mat, pair_exprsn1, pair_exprsn2, all_cells):
#     print(lr_index)

    p1_mean = pair_exprsn1.iloc[:,lr_index].mean()
    p2_mean = pair_exprsn2.iloc[:,lr_index].mean()
    dfa = pair_exprsn1[pair_exprsn1.iloc[:,lr_index] >p1_mean].iloc[:,lr_index]
    dfb = pair_exprsn2[pair_exprsn2.iloc[:,lr_index] >p2_mean].iloc[:,lr_index]
    thresh = round(pair_exprsn2.shape[0]*0.05)
    if len(dfa) > thresh:
        dfa = dfa.astype(float).nlargest(thresh)
    if len(dfb) > thresh:
        dfb = dfb.astype(float).nlargest(thresh)
    if len(dfa) < len(dfb):
        adj_mat = pair_sum_over_dfa(adj_mat, dfa, dfb, all_cells)
    else:
        adj_mat = pair_sum_over_dfb(adj_mat, dfa, dfb, all_cells)
    return adj_mat

#     p1_mean = pair_exprsn1.iloc[:,lr_index].mean()
#     p2_mean = pair_exprsn2.iloc[:,lr_index].mean()
#     dfa = pair_exprsn1[pair_exprsn1.iloc[:,lr_index] >p1_mean].iloc[:,lr_index]
 #     dfb = pair_exprsn2[pair_exprsn2.iloc[:,lr_index] >p2_mean].iloc[:,lr_index]
    #     # thresh = round(len(pair_exprsn1.iloc[:,lr_index])*0.1)
# def gene_pair_score(lr_index, adj_mat):
# #     print(lr_index)
#     p1_mean = pair_exprsn1.iloc[:,lr_index].mean()
#     p2_mean = pair_exprsn2.iloc[:,lr_index].mean()
#     dfa = pair_exprsn1[pair_exprsn1.iloc[:,lr_index] >p1_mean].iloc[:,lr_index]
#     dfb = pair_exprsn2[pair_exprsn2.iloc[:,lr_index] >p2_mean].iloc[:,lr_index]
#     thresh = 1000
#     if len(dfa) > thresh:
#         dfa = dfa.astype(float).nlargest(thresh)
#     if len(dfb) > thresh:
#         dfb = dfb.astype(float).nlargest(thresh)
#     if len(dfa) < len(dfb):
#         adj_mat = pair_sum_over_dfa(adj_mat, dfa, dfb)
#     else:
#         adj_mat = pair_sum_over_dfb(adj_mat, dfa, dfb)
#     return adj_mat
# sum_pair_scores = np.sum(all_pair_scores, axis = 0)




