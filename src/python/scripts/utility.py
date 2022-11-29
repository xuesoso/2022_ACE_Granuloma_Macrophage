# Path: src/scripts/utility.py
# Common utility functions used by analysis scripts
# Author: Yuan Xue (xuesoso@gmail.com)
# Github: https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
# Updated on 28 Nov 2022
# License: MIT

#### Import variables and libraries ####
from ..variables.granuloma_lib import *

def save_h5ad(ad, filename, compression='gzip'):
    """
    Save AnnData object to h5ad file

    Parameters
    ----------
    ad : AnnData
        AnnData object
    filename : str
        Path to output file
    compression : str
        Compression method
    """
    #### Write file atomically to avoid incomplete file ####
    temp_filename = ('/'.join(filename.split('/')[:-1]) + '/__temp__'
                        + filename.split('/')[-1])
    ad.write_h5ad(temp_filename, compression=compression)
    os.rename(temp_filename, filename)

def transfer_ad_key(obs_to, obs_from, replace=False):
    """
    Transfer the column values from one pd.DataFrame object to another

    Parameters
    ----------
    obs_to : pd.DataFrame
        The pd.DataFrame object to transfer the column values to
    obs_from : pd.DataFrame
        The pd.DataFrame object to transfer the column values from
    replace : bool
        Whether to replace the existing column values in obs_to
    """
    for key in obs_from.keys():
        if key not in obs_to.keys():
            obs_to[key] = obs_from[key].copy()
        elif replace and key in obs_to.keys():
            obs_to[key] = obs_from[key].copy()

def load_marker_gene(directory):
    """
    Load marker genes from a directory. Presumes that each file is named with an
    underscore as the delimiter. The first part is the cluster name.

    Parameters
    ----------
    directory : str
        Path to marker gene directory

    Returns
    -------
    marker_gene : dict
        A dictionary of marker genes
    """
    fns = glob(directory+'/*')
    marker = {}
    for fn in fns:
        cl = fn.split('/')[-1].split('_')[0]
        marker[cl] = pd.read_csv(fn, sep='\t', index_col=0)
    return marker

def load_cell_cycle_genes(fn=infd+'../external/mouse_cell_cycle_marker_genes.csv'):
    """
    Load cell cycle marker genes

    Parameters
    ----------
    fn : str
        Path to cell cycle marker gene file

    Returns
    -------
    s_genes : list
        List of S-phase marker genes
    g2m_genes : list
        List of G2M-phase marker genes
    """
    ccg = pd.read_csv(fn)
    s_genes = [x.strip().capitalize() for x in ccg.loc[[True if 'S' in x else
        False for x in ccg['Phase']]]['Gene'].values]
    g2m_genes = [x.strip().capitalize() for x in ccg.loc[[True if 'G2' in x or
        'M' in x else False for x in ccg['Phase']]]['Gene'].values]
    return s_genes, g2m_genes

def load_panglaoDB(ad, processed_panglao=None):
    """
    Load PanglaoDB cell type annotation

    Parameters
    ----------
    ad : AnnData
        AnnData object
    processed_panglao : str
        Path to processed PanglaoDB file

    Returns
    -------
    panglaoDB_marker : pd.DataFrame
        A dataframe of panglaoDB marker genes and their frequency in each cell
        type
    """
    if processed_panglao is None:
        processed_panglao = infd+'../external/processed_PanglaoDB_markers_24_Jan_2020.tsv.gz'
    if os.path.isfile(processed_panglao):
        print('reading existing {:}'.format(processed_panglao))
        panglaoDB_marker = pd.read_csv(processed_panglao, sep='\t')
    else:
        panglaoDB_marker = pd.read_csv(infd+'../external/PanglaoDB_markers_24_Jan_2020.tsv.gz', sep='\t')
        panglaoDB_marker['GeneSymbol'] = [x.capitalize() for x in panglaoDB_marker['official gene symbol']]
        ind = np.in1d(ad.var_names, panglaoDB_marker['GeneSymbol'])
        e = ad[:, ind].X.sum(axis=0).A.flatten() > 0
        ind = np.in1d(panglaoDB_marker['GeneSymbol'], ad.var_names[ind][e])
        panglaoDB_marker = panglaoDB_marker.loc[ind]
        panglaoDB_marker['frequency'] = [(panglaoDB_marker['GeneSymbol'] == x).sum() for x in panglaoDB_marker['GeneSymbol']]
        panglaoDB_marker.to_csv(processed_panglao, sep='\t', compression='gzip')
    return panglaoDB_marker

def load_GSOA(outDir):
    """
    Load GSOA results

    Parameters
    ----------
    outDir : str
        Path to GSOA output directory

    Returns
    -------
    gsoa : pd.DataFrame
        GSOA results
    """
    out = {}
    fns = glob(os.path.join(outDir, '*.tsv.gz'))
    for fn in fns:
        lastl = fn.split('/')[-1]
        group = lastl.split('_')[0]
        if group not in out:
            out[group] = {}
        geneset = '_'.join(lastl.split('_')[1:]).split('.tsv.gz')[0]
        out[group][geneset] = pd.read_csv(fn, sep='\t', index_col=0)
    return out

def read_msigdb(fn, capitalize=True):
    """
    Load MSigDB gene sets from GSEA

    Parameters
    ----------
    fn : str
        Path to MSigDB file
    capitalize : bool
        Whether to capitalize gene names

    Returns
    -------
    out : dict
        A dictionary of MSigDB gene sets
    """
    capitalize=True
    out = {}
    with open(fn) as f:
        for l in f:
            spl = l.strip().split()
            members = spl[3:]
            if capitalize:
                members = [x.capitalize() for x in members]
            out[spl[0]] = members
    return out

def get_combined_goa(X, excluded_pathways=[], upregulated=True, intersect_threshold=3, pthreshold=0.1):
    """
    Retrieve gene set over-representation analysis (GSOA) results

    Parameters
    ----------
    X : AnnData
        AnnData object
    excluded_pathways : list
        List of pathways to exclude
    upregulated : bool
        Whether to use upregulated genes
    intersect_threshold : int
        Minimum number of genes in intersection
    pthreshold : float
        p-value threshold

    Returns
    -------
    combined : pd.DataFrame
        Combined GOA results
    """

    combined = pd.DataFrame()
    for k in X:
        if upregulated:
            if k not in excluded_pathways and k[-2:] == 'Up':
                combined = pd.concat([X[k], combined])
        else:
            if k not in excluded_pathways and k[-4:] == 'Down':
                combined = pd.concat([X[k], combined])
    combined = combined[(combined['intersect_size'] >= intersect_threshold)]
    if len(combined) > 0:
        combined['combined_fdr'] = sm.stats.multipletests(combined['pval'], method='fdr_bh')[1]
        combined = combined.sort_values('combined_fdr')
        combined = combined[combined['combined_fdr'] < pthreshold]
        combined['-log10_fdr'] = -np.log10(combined['combined_fdr'])
    return combined

def from_adata_to_df(ad, obs=None, var=None):
    """
    Convert anndata object to a dataframe

    Parameters
    ----------
    ad : anndata object
        anndata object to convert
    obs : list of str, optional
        list of observation names to include in the dataframe, by default None
    var : list of str, optional
        list of variable names to include in the dataframe, by default None

    Returns
    -------
    df : pd.DataFrame
        dataframe with the data
    o1 : additional dataframe.
    o2 : additional dataframe.
    """
    df = pd.DataFrame(ad.X.A, index=ad.obs_names, columns=ad.var_names)
    if obs is not None:
        o1 = ad.obs[obs]
    else:
        o1 = None
    if var is not None:
        o2 = ad.var[var]
    else:
        o2 = None
    return df, o1, o2

def calculate_vector_mean(V, ignore_zeros=False, axis=None, operation='mean'):
    """
    Calculate the vector mean of a matrix

    Parameters
    ----------
    V : np.array
        matrix to calculate the vector mean
    ignore_zeros : bool, optional
        whether to ignore zeros, by default False
    axis : int, optional
        axis to calculate the mean, by default None
    operation : str, optional
        operation to perform, by default 'mean'. Options include "median" or
        "mean"

    Returns
    -------
    t : np.array
        vector mean
    """
    if ignore_zeros == False:
        if operation == 'mean':
            return np.mean(V, axis=axis)
        elif operation == 'median':
            return np.median(V, axis=axis)
        else:
            raise ValueError('{:} is not a valid operation'.format(operation))
    else:
        t = V.copy()
        t[t <= 0] = np.nan
        if operation == 'mean':
            t = np.nanmean(t, axis=axis)
        elif operation == 'median':
            t = np.nanmedian(t, axis=axis)
        else:
            raise ValueError('{:} is not a valid operation'.format(operation))
        t = np.nan_to_num(t)
        return t

def average_by_clusters(mat,
                        clusters,
                        return_clusters=False,
                        ignore_zeros=False,
                        operation='mean'
                        ):
    """
    Function to convert an expression matrix (n x m) to a expression-average
    matrix by a list of clusters with a set size, c.


    Parameters
    ----------
    mat : pd.DataFrame
        Expression matrix (n x m)
    clusters : pd.Series
        Series of clusters (n)
    return_clusters : bool, optional
        Return a matrix in the shape of (c x m) if True, by default False
    ignore_zeros : bool, optional
        Calculate non-zero average if True, by default False
    operation : str, optional
        Operation to use when calculating the average, by default 'mean'

    Returns
    -------
    new_mat : pd.DataFrame
        Expression average matrix
    """
    clusters = np.array(clusters)
    unique, counts = np.unique(clusters, return_counts=True)
    return_df = False
    if isinstance(mat, pd.DataFrame):
        col_names = mat.columns.values
        row_names = mat.index.values
        mat = mat.T.values
        return_df = True
    else:
        col_names = np.arange(mat.shape[0])
        row_names = np.arange(mat.shape[1])
    if len(np.shape(mat)) < 2:
        if return_clusters == False:
            new_mat = np.zeros(np.shape(clusters)[0])
            for i in unique:
                new_mat[clusters == i] = calculate_vector_mean(
                    mat[clusters == i], ignore_zeros=ignore_zeros,
                    operation=operation)
        else:
            new_mat = np.zeros(np.shape(unique)[0])
            for count, i in enumerate(unique):
                new_mat[count] = calculate_vector_mean(
                    mat[clusters == i], ignore_zeros=ignore_zeros,
                    operation=operation)

    else:
        if return_clusters == False:
            new_mat = np.zeros(np.shape(mat))
            size_of_mat = new_mat.shape[0]
            for i, cnt in zip(unique, counts):
                new_mat[:, clusters == i] = np.repeat(
                    calculate_vector_mean(mat[:, clusters == i],
                                          ignore_zeros=ignore_zeros,
                                          axis=1, operation=operation),
                    cnt).reshape(size_of_mat, cnt)
        else:
            new_mat = np.zeros((np.shape(mat)[0], np.shape(unique)[0]))
            count = 0
            for i, cnt in zip(unique, counts):
                new_mat[:, count] = calculate_vector_mean(
                    mat[:, clusters == i], axis=1, ignore_zeros=ignore_zeros,
                    operation=operation)
                count += 1
    if return_df == True or return_clusters == True:
        if return_clusters == True:
            new_mat = pd.DataFrame(new_mat.T, index=unique, columns=col_names)
        else:
            new_mat = pd.DataFrame(new_mat.T, \
                                   index=row_names, columns=col_names)
    return new_mat
