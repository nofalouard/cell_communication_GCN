import pandas as pd
from kneed import KneeLocator
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import scipy.stats as stats
import mxnet.ndarray as nd
from itertools import compress
from models import build_features, build_model
from layers import SpectralRule
import mxnet as mx


def find_clusters(adata, df):
    # inds=df.columns[0]
    # df.set_index(inds, inplace=True)  
    # cell_rows = list(df.index)
    # cell_cols = list(df.columns)
    # uni_rows = list(set(cell_rows) - set(cell_cols))
    # uni_cols = list(set(cell_cols) - set(cell_rows)) 
    # new_df = pd.DataFrame(0.0, index = uni_cols, columns = uni_rows)
    # adj_df = df.append(new_df)
    # adj_df = adj_df.fillna(0)
    # adj_df = adj_df.reindex(columns = adj_df.index)
    # adj = adj_df.to_numpy()
    
    # inds=df.columns[0]
    # df.set_index(inds, inplace=True)
    A = df.to_numpy()
    print(A.shape)
    A = stats.zscore(A, axis=None, ddof=1)

    # in expression matrix, keep only the cells selected in the gcn adjacency matrix 
    all_cells = list(adata.obs.index)
    filt_cells = list(df.index)

    print("there are " + str(len(all_cells)) + " total cells and " + str(len(filt_cells)) + " filtered cells")
    index_dict = dict((value, idx) for idx,value in enumerate(all_cells))
    to_filt = [index_dict[x] for x in filt_cells]
    exprsn= adata.X.A
    
    exprsn_gcn = exprsn[to_filt,:]
    print("the filtered exprsn matrix has dims " + str(exprsn_gcn.shape))
    pca = PCA(n_components=30)
    pca.fit(exprsn_gcn)
    embed=pca.transform(exprsn_gcn)
    X_pca = embed
    print("building model using pca as input features")
    # build gcn model using the adjacency matrix and input features
    model_1, features_1 = build_model(nd.array(A), nd.array(X_pca))
    feats = model_1(nd.array(X_pca))
    feats = feats.asnumpy()
    print("beginning to find optiml number of clusters")
    # find the optimal number of clusters using elbow method
    # A list holds the SSE values for each k
    sse = []
    for k in range(1, 11):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(feats)
        sse.append(kmeans.inertia_)
    
    kl = KneeLocator(
    range(1, 11), sse, curve="convex", direction="decreasing")

    opt_clusts = kl.elbow
    print("the optimal number of clusters is ")
    print(opt_clusts)

    # apply the kmeans algorithm to cluster the outputted feature matrix and find the optimal number of clusters for the 
    # dataset
    kmeans = KMeans(n_clusters=opt_clusts)
    # feats = feats.asnumpy()
    kmeans.fit(feats)
    all_predictions = kmeans.predict(feats)

    # identify which cells belong to each cluster using the index of each predicted label
    # store the cell labels that belong to each cluster in a dataframe

    clust_dict = {}
    for i in range(0, len(range(opt_clusts))):
        clust_i = list(compress(filt_cells,all_predictions==i ))
        print("the length of this cluster is")
        print(len(clust_i))
        clust_dict["Cluster_" + str(i)]= clust_i

    clust_df = pd.DataFrame({ key:pd.Series(value) for key, value in clust_dict.items() })
    cell_annot = adata.obs.loc[filt_cells]
    exprsn_filt = pd.DataFrame(exprsn_gcn, index=filt_cells, columns= adata.var_names)
    
    return clust_df, cell_annot, exprsn_filt