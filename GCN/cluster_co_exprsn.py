import itertools
import random
import numpy as np
import pandas as pd
from itertools import compress
from itertools import combinations
from joblib import Parallel, delayed
from scipy import stats
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt

def find_pair_coexprsd(clust_df, clust_annot, pair, adata):
    clust_inds = clust_annot.columns[0]
    clust_annot.set_index(clust_inds, inplace=True)
    
    clust_df = clust_df.drop(columns=clust_df.columns[0])
    exprsn = pd.DataFrame(adata.X.A)
    exprsn.set_index(adata.obs.index, inplace=True)
    exprsn.columns = adata.var_names

    clust_cells = {}
    for (name, data) in clust_df.iteritems():
        clust_cells[name] = data.dropna()
    
    clust_labels = {}
    for key, item in clust_cells.items():
        clust_labels[key] = np.array(clust_annot.loc[item, 'assign.level3_anno'], dtype = 'str')
    
    ctype_combs = {}
    for key, item in clust_labels.items():
        types_lst = list(set(item))  
        ctype_combs[key] = list(combinations(types_lst, 2)) + [(x,x) for x in types_lst]
    # calculate pearson correlation coefficients for all gene pairs within each cluster 

    all_coexprsd = []
    intersect = []
    print("here")
    all_coexprsd = []
    intersect = []
    pair_found = False
    for key, clust in clust_cells.items():
        ctypes = list(set(clust_labels[key]))
        label_dict = {ctypes[i]: list(compress(clust,clust_labels[key]== ctypes[i])) for i in range(len(ctypes))}
        if pair[0] in list(label_dict.keys()) and pair[0] in list(label_dict.keys()):
            pair_found = True
            #for pair in ctype_combs[key][0:1]:
            pair_cells = label_dict[pair[0]] + label_dict[pair[1]]
            corr_array = pair_corr_analysis(exprsn, label_dict, pair)
            print('got corr array')
            sig_pairs_names = clust_coexprsd(exprsn,pair_cells, corr_array)
            intersect = list(set(all_coexprsd) & set(sig_pairs_names))
            all_coexprsd = list(set(sig_pairs_names + all_coexprsd))
        else:
            print(key + " doesnt have both cell types")

    print("union of pair co-expressed genes: " + str(len(all_coexprsd)))
    print("intersect of pair co-expressed genes: " + str(len(intersect)))

    return all_coexprsd, pair_found

def find_clust_coexprsd(clust_df, adata):
    clust_df = clust_df.drop(columns=clust_df.columns[0])
    exprsn = pd.DataFrame(adata.X.A)
    exprsn.set_index(adata.obs.index, inplace=True)
    exprsn.columns = adata.var_names

    clust_cells = {}
    for (name, data) in clust_df.iteritems():
        clust_cells[name] = data.dropna()
    
    # calculate pearson correlation coefficients for all gene pairs within each cluster 

    all_coexprsd = []
    intersect = []
    print("here")
    for key, clust in clust_cells.items():
        clust_corr = clust_corr_analysis(exprsn, clust)
        print("finished with obs corrs")
        sig_pairs_names = clust_coexprsd(exprsn,clust, clust_corr)
        # store the significantly co-expressed gene pairs in the all_coexprsd list
        intersect = list(set(all_coexprsd) & set(sig_pairs_names))
        all_coexprsd = list(set(sig_pairs_names + all_coexprsd))
        # print(len(all_coexprsd), len(intersect))
    print("union of cluster co-expressed genes: " + str(len(all_coexprsd)))
    print("intersect of cluster co-expressed genes: " + str(len(intersect)))

    return all_coexprsd
    
'''
   Identifies most significantly co-expressed gene pairs in a cell cluster. 

   Input:
   - Gene expression matrix for whole dataset
   - Numpy array with correlation coefficients for each gene pair's expression within that cluster. 

   Conducts permutation test on gene labels to generate a null distribution of correlation coefficients. Uses the 
   distribution to then calculate the p value of each gene pair's observed correlation. 
'''
def clust_coexprsd(exprsn, clust, clust_corr):
    num_nan = np.count_nonzero(np.isnan(clust_corr))
    corr_sort = np.argsort(clust_corr, axis = None)
    top_corrs = corr_sort[round(0.95*(len(corr_sort)-num_nan)):-num_nan]
    print("there are " + str(len(top_corrs)) + " top correlated pairs")
    num_genes = exprsn.shape[1]
    top_rows = get_rows(top_corrs, num_genes)
    top_cols = get_cols(top_corrs, num_genes)
    obs_corrs=clust_corr[top_rows,top_cols]
    dist = permuted_genes_corr(exprsn, clust, 2, [obs_corrs])
    print("finished null distribution analysis")
    zscores = stats.zscore(dist, axis = 1, nan_policy='omit' )
    pos_scores = zscores[0] > 0
    sig_scores = stats.norm.sf(abs(zscores[0])) < 0.05
    sig_inds = np.where(pos_scores & sig_scores)[0]
    sig_pairs_flat = top_corrs[sig_inds] 
    sig_pairs_rows = get_rows(sig_pairs_flat, num_genes) 
    sig_pairs_cols = get_cols(sig_pairs_flat, num_genes) 
    sig_pairs = list(zip(sig_pairs_rows, sig_pairs_cols))
    gene_lst = list(exprsn.columns)
    sig_pairs_names = list(map(lambda x: (gene_lst[x[0]],gene_lst[x[1]]), sig_pairs))
    print("there are " + str(len(sig_pairs_names)) + " sig pairs in this clust")
    return sig_pairs_names

'''
   Permutes the gene labels of an expression matrix n_iter times and calculates the pearson correlation coefficients
   for those genes each time. Returns a 2D array with the resulting correlation coefficients 
   between genes whose index belongs in the gene_rows list and genes whose index belongs in gene_cols.

'''

def permuted_genes_corr(exprsn, clust, n_iter, obs_corrs):
    for i in range(n_iter):
        cells_df = exprsn.loc[clust]
        cells_array = np.array(cells_df)
        np.apply_along_axis(np.random.shuffle,1,cells_array)
        corr_array = np.corrcoef(cells_array, rowvar=False)
        # select 25 random values from corr_array
        sample_size = 25
        rand_rows = np.random.choice(corr_array.shape[0], sample_size,replace=False)
        rand_cols = np.random.choice(corr_array.shape[0], sample_size,replace=False)
        rand_corrs = corr_array[rand_rows,rand_cols]
        repeat_corrs = np.tile(rand_corrs, (obs_corrs[0].shape[0],1))
        repeat_corrs = np.transpose(repeat_corrs)
        obs_corrs = np.append(obs_corrs,repeat_corrs, axis = 0)
    return obs_corrs


'''
   Convert a list of 1D array indices to the row indices of its reshaped, 2D array form, with number of columns
   defined by num_cols
'''
def get_rows(flat_ind, num_cols):
    return np.trunc(flat_ind/num_cols).astype('int')

'''
   Convert a list of 1D array indices to the column indices of its reshaped, 2D array form, with number of columns
   defined by num_cols
'''
def get_cols(flat_ind, num_cols):
    return flat_ind%num_cols


'''
   Input: Expression matrix and list of selected cells
   Does: Based on selected cells' gene expression values, calculates pearson correlation coefficient 
         between each gene pair
   Returns: pearson correlation matrix
'''
def clust_corr_analysis(exprsn, clust):
    cells_df = exprsn.loc[clust]
    cells_array = np.array(cells_df)
    corr_array = np.corrcoef(cells_array, rowvar=False)
    corr_array[list(range(corr_array.shape[0])), list(range(corr_array.shape[0]))] = np.nan
    return corr_array


'''
Input: Expression matrix, dictionary listing which cells in matrix belong to each cell type, and a selected 
pair of cell types
Does: computes pearson correlation matrix for cells with a type that is specified in the "pair" argument
Returns: pearson correlation matrix
'''
def pair_corr_analysis(exprsn, label_dict, pair):
    cells = label_dict[pair[0]] + label_dict[pair[1]]
    cells_df = exprsn.loc[cells]
    cells_array = np.array(cells_df)
    corr_array = np.corrcoef(cells_array, rowvar=False)
    corr_array[list(range(corr_array.shape[0])), list(range(corr_array.shape[0]))] = np.nan
    return corr_array