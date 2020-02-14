# Imports
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import issparse
import numpy as np
from anndata import AnnData

from .utils import _calc_DoTscore, _sim_DoTscore, _randomde, _mean_var, _calc_DoTscore_matrix

# User functions
def custom_scale(adata, mean = None):
    """Function to scale the expression data by mean and variance - (counts-mean(counts))/sd(counts). Based on the scanpy function for compatibility.

    Parameters
    ----------
    adata : AnnData
        The AnnData object which .X slot is used for calculation
    mean : array
        A numpy array with mean gene expression values to be used for scaling. If None, global mean will be used.

    Returns
    -------
    Updates the expression values in place

    """
    if isinstance(adata, AnnData):
        if issparse(adata.X):
            adata.X = adata.X.toarray()
        custom_scale(adata.X, mean = mean)
        return None

    X = adata

    meanTOTAL, var = _mean_var(X)
    scale = np.sqrt(var)

    if mean is None:
        mean = meanTOTAL

    X -= mean
    scale[scale == 0] = 1e-12
    X /= scale

def get_DoTscore(adata,
                 de,
                 allgenes = None,
                 allfolds = None,
                 simno = 500,
                 id_col = 'target',
                 weight_col = 'log2FoldChange',
                 zscore = False):
    """ Generated a z-score version of the FCscore. For each cells the scaled expression values
    are multiplied by fold change supplied in de (only if present in that data frame)
    and summed up. This process is also simmulated with random genes and random fold
    changes and a z-score is calculated.

    Parameters
    ----------
    adata : AnnData
        The AnnData object which .X slot is used for calculation
    de : PandasDataFrame
        DataFrame with differentially expressed genes, default columns used are: 'target' and 'log2FoldChange'
    allgenes : numpy array or pandasSeries
        numpy array with all the genes expressed in the cells where de is generated
    allfolds : numpy array or pandasSeries
        an array of fold changes from which the simulation should draw
    simno : int
        number of simulations to run (for estimation of z score)
    id_col : str
        name of the column in de which contains the gene ids
    weight_col : str
        name of the column in de which contains the weights (default: log2FoldChange)

    Returns
    -------
    type pandasSeries
        z-scores of FCscores

    """

    if isinstance(adata, AnnData): pass
    else: raise ValueError('Data needs to be an Anndata object')


    size = len(de.index)

    score = _calc_DoTscore(adata,
                           de,
                           id_col = id_col,
                           weight_col = weight_col)

    if zscore:

        if isinstance(allgenes, np.ndarray): pass
        else: allgenes = allgenes.values
        if isinstance(allfolds, np.ndarray): pass
        else: allfolds = allfolds.values

        stat = _sim_DoTscore(adata,
                            allgenes,
                            allfolds,
                            size = size,
                            simno = simno)
        
        score = (score - stat['mean']) / stat['sd']

    return(score)


def qfilt(values, qmax = 0.99, qmin = 0.01):
    """Quantile value filter. Changes all values above a certain quantile
    to the value equal to that quantile

    Parameters
    ----------
    values : PandasSeries, numpy array, list or tuple
        List of values to be filtered
    qmax : numeric
        The upper quantile used for filtering (all values above will be equal
        to the value of this quantile)
    qmin : numeric
        The lower quantile used for filtering (all values below will be equal
        to the value of this quantile)

    Returns
    -------
    type : PandasSeries
        Series with filtered values

    """
    if isinstance(values, pd.Series): pass
    else: values = pd.Series(values)

    values2 = values.copy()
    vmax = values2.quantile(q=qmax)
    vmin = values2.quantile(q=qmin)

    values2.loc[values2 > vmax] = vmax
    print("Filtering: " + str(len(values2.loc[values > vmax])) + ' cells above threshold' + ' (' + str(vmax) + ')')
    values2.loc[values2 < vmin] = vmin
    print("Filtering: " + str(len(values2.loc[values < vmin])) + ' cells below threshold' + ' (' + str(vmin) + ')')

    return(values2)

def get_genescore_pergroup(adata,
                          de,
                          id_col = 'target',
                          weight_col = 'log2FoldChange',
                          group = 'leiden',
                          sortby = '0',
                          gene_symbols = None):
    """ Calculates DoT score contribution coming from each gene, averaged for each cell group (e.g. clusters)

    Parameters
    ----------
    adata : AnnData
        AnnData object with the expression values in the .X slot and the group
        column in the .obs slot
    de : PandasDataFrame
        DataFrame with differentially expressed genes, default columns used are: 'target' and 'log2FoldChange'
    id_col : str
        Name of the column in de with gene ids
    id_col : str
        name of the column in de which contains the gene ids
    weight_col : str
        name of the column in de which contains the weights (default: log2FoldChange)
    group : str
        name of the column in the .obs slot of AnnData object with cell groups (needs to contain categorical data)
    sortby : str
        name of the cell group by which the genes should be sorted (based on scores)
    gene_symbols : str
        If not none, indicates the name of the column in the .var slot of the
        AnnData object which contains gene symbols

    Returns
    -------
    type : pandasDataFrame
        DataFrame with FCscores for each group, summed per cell group

    """

    df = pd.DataFrame(index = adata.var.index)

    for i in adata.obs[group].values.categories:

        sub = _calc_DoTscore_matrix(adata = adata[adata.obs[group] == i,:], de = de, id_col = id_col, weight_col = weight_col)
        # sub = x[adata.obs[group] == i,:]
        sub = sub.mean(axis = 0)
        df[i] = sub

    df.columns = df.columns.astype('object')
    
    if not(gene_symbols is None):
        df['symbol'] = adata.var[gene_symbols]

    df = df.sort_values(by = sortby, axis = 0, ascending = True)
    # df = df.sort_index(axis = 1)

    return(df)

def cmap_RdBu(values, vmin = None, vmax = None):
    """Generates a blue/red colorscale with white value centered around the value 0

    Parameters
    ----------
    values : PandasSeries, numpy array, list or tuple
        List of values to be used for creating the color map
    vmin : type
        Minimum value in the color map, if None then the min(values) is used
    vmax : type
        Maximum value in the color map, if None then the max(values) is used

    Returns
    -------
    type
        Description of returned object.

    """
    if vmin != None:
        scoremin = vmin
    else:
        scoremin = min(values)
    if vmax != None:
        scoremax = vmax
    else:
        scoremax = max(values)

    from matplotlib.colors import LinearSegmentedColormap
    cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, 'blue'),
                                                        (-scoremin/(scoremax-scoremin), 'white'),
                                                        (1, 'red')])
    return(cmap2)
