import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import issparse
import numpy as np
from anndata import AnnData

def _calc_DoTscore(adata, de, id_col = 'target', weight_col = 'log2FoldChange'):
    """ Calculates DoT scores per cell, summing up all the FCscores for each gene
    that is differentially expressed.

    Parameters
    ----------
    adata : AnnData
        The AnnData object which .X slot is used for calculation
    de : PandasDataFrame
        DataFrame with differentially expressed genes, default columns used are: 'target' and 'log2FoldChange'
        log2FoldChange column
    id_col : str
        name of the column in de which contains the gene ids
    weight_col : str
        name of the column in de which contains the weights (default: log2FoldChange)

    Returns
    -------
    type
        Description of returned object.

    """

    #Generating a vector with weights and 0 for not changing genes
    df = pd.DataFrame({'id' : adata.var.index.tolist()})
    df = df.merge(de[[id_col, weight_col]], left_on='id', right_on = id_col, how= 'left')
    NAs = np.isnan(df[weight_col]).tolist()
    df.loc[NAs, weight_col] = 0

    scores = df[weight_col].values

    if issparse(adata.X):
        X = adata.X.toarray()
    else: X = adata.X

    cellscores = np.matmul(X, scores)
    return(cellscores)

def _sim_DoTscore(adata, allgenes, allfolds, size, simno):
    """Simulates DoT score given a set of genes and weights

    Parameters
    ----------
    adata : AnnData
        The AnnData object which .X slot is used for calculation
    de : PandasDataFrame
        DataFrame with differentially expressed genes, default columns used are: 'target' and 'log2FoldChange'
        log2FoldChange column
    allgenes : numpy array
        numpy array with all the genes expressed in the cells where de is generated
    allfolds : numpy array
        an array of fold changes from which the simulation should draw
    size : int
        number of non-zero weights (typically number of DE genes)
    simno : int
        number of simulations to run (for estimation of z score)

    Returns
    -------
    type : PandasDataFrame
        DataFrame with the mean and standard deviation statistic per cell

    """
    score = np.zeros((len(adata.obs.index), simno))
    folds = np.zeros((len(adata.var.index), simno))
    folds = pd.DataFrame(folds, index = adata.var.index)

    for i in range(0, simno):
        de1 = _randomde(allgenes, allfolds, size = size)
        #Subsetting only for the genes present in the adata dataset
        sub = de1['id'].isin(adata.var.index)
        de1 = de1.loc[sub]

        folds.loc[de1['id'].tolist(), i] = de1['weights'].tolist()

    folds = np.array(folds)

    if issparse(adata.X):
        X = adata.X.toarray()
    else: X = adata.X

    cellscores = np.matmul(X, folds)

    meanscore = np.mean(cellscores, 1)
    sdscore = np.var(cellscores, 1)**0.5
    statdf = pd.DataFrame({'mean' : meanscore, 'sd' : sdscore}, index = adata.obs.index)
    return(statdf)

def _randomde(allgenes,
              allfolds,
              size):
    """Randomly select genes from the allgenes array and fold changes from the
    allfolds array. Size argument indicates how many to draw.

    Parameters
    ----------
    allgenes : numpy array
        numpy array with all the genes expressed in the cells where de is generated
    allfolds : numpy array
        an array of fold changes from which the simulation should draw
    size : int
        number of non-zero weights (typically number of DE genes)

    Returns
    -------
    type : PandasDataFrame
        DataFrame with randomly chosen genes and weights.

    """
    rdgenes = np.random.choice(allgenes, size, replace = False)
    rdfolds = np.random.choice(allfolds, size, replace = False)
    rdDF = pd.DataFrame({'id' : rdgenes, 'weights' : rdfolds})
    return(rdDF)

def _mean_var(counts):
    #Adapted from scanpy to keep numbers exactly the same
    mean = counts.mean(axis=0)
    
    if issparse(counts):
        sq = counts.multiply(counts)
        #If using sparse matrices the need to be flattened (A1)
        mean_sq = sq.mean(axis = 0).A1
        mean = mean.A1
    else:
        sq = np.multiply(counts, counts)
        mean_sq = sq.mean(axis = 0)

    n = counts.shape[0]
    var = (mean_sq - mean**2) * (n/(n-1))
    return mean, var

def _calc_DoTscore_matrix(adata,
                          de,
                          id_col = 'target',
                          weight_col = 'log2FoldChange'):
    """ Generates a cells x genes FCscore matrix, which contains multiplied values
    if scaled expression values and the observed fold change. 

    Parameters
    ----------
    adata : AnnData
        the AnnData object used as a source of expression values. .X slot is used.
    de : PandasDataFrame
        Data.frame with a column with gene ids (name set in id_col) and fold changes in a column called log2FoldChange (default)
    id_col : str
        name of the column in the de to used as gene ids
    weight_col : str
        name of the column in de with fold change values

    Returns
    -------
    type : numpy array
        Array cells x genes with scores

    """

    df = pd.DataFrame(data = adata.var.index.tolist())
    df = df.merge(de[[id_col, weight_col]], left_on=0, right_on = id_col, how= 'left')

    NAs = np.isnan(df[weight_col]).tolist()
    df.loc[NAs, weight_col] = 0

    des = df[weight_col].values

    scores = np.zeros(adata.shape)
    #' Start with all 0s and update positions just for the genes with non0 values in des (should be faster)

    toupdate = np.where(df[weight_col] != 0)

    for i in toupdate:
        if(adata.X.ndim > 1): scores[:,i] = adata.X[:,i] * des[i]
        else: scores[:,i] = adata.X[i] * des[i]

    return(scores)
