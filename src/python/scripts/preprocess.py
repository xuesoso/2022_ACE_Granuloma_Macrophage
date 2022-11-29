# Path: src/scripts/preprocess.py
# Define functions used for data normalization
# Author: Yuan Xue (xuesoso@gmail.com)
# Github: https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
# Updated on 28 Nov 2022
# License: MIT

#### Import variables and libraries ####
from ..variables.granuloma_lib import *


def normalize_per_cell(adata, sum_norm='cell_median', max_fraction_genes=0):
    """
    Normalize scRNA-seq data by library sequencing depth.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    sum_norm : `str` (default: 'cell_median')
        Normalization method for total counts per cell.
        If `cell_median`, normalize each cell by the sum of its read counts
        and multiply by the median of read count sums across all samples.
        If `cell_mean`, normalize each cell by the sum of its read counts
        and multiply by the mean of read count sums across all samples.
    max_fraction_genes : `float` (default: 0)
        If `max_fraction_genes` is not 0, then only genes with a fraction of
        total counts greater than or equal to `max_fraction_genes` are
        considered for the normalization.
    """
    D = adata.X
    if isinstance(D, np.ndarray):
        D=sp.sparse.csr_matrix(D,dtype='float32')
    else:
        if str(D.dtype) != 'float32':
            D=D.astype('float32')
        D.sort_indices()

    if(D.getformat() == 'csc'):
        D=D.tocsr();

    if max_fraction_genes > 0:
        fraction = D / D.sum(1)
        mask = (fraction) < 0.05

    if (sum_norm == 'cell_median'):
        if max_fraction_genes > 0:
            s = D[mask].sum(1).A.flatten()
        else:
            s = D.sum(1).A.flatten()
        sum_norm = np.median(s)
        D = D.multiply(1 / s[:,None] * sum_norm).tocsr()

    elif (sum_norm == 'cell_mean'):
        if max_fraction_genes > 0:
            s = D[mask].sum(1).A.flatten()
        else:
            s = D.sum(1).A.flatten()
        sum_norm = np.mean(s)
        D = D.multiply(1 / s[:,None] * sum_norm).tocsr()
    adata.X = D

def normalize(adata, norm='log2', div=1, **args):
    """
    Transform the depth-normalized read count matrix to different space.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    norm : `str` (default: 'log2')
        Normalization method.
        If `log2`, (log2+1) transform the normalized read count matrix.
        If `log1p`, (log+1) transform the normalized read count matrix.
        If `pearson`, Pearson residuals of the normalized read count matrix.
    args : `dict`
        Arguments for :func:`pearson_residuals`.
    """
    norm = norm.lower()
    D = adata.X

    if norm is None:
        D.data[:] = (D.data / div)

    elif(norm == 'log1p'):
        D.data[:] = np.log(D.data / div + 1)

    elif(norm == 'log2'):
        D.data[:] = np.log2(D.data / div + 1)

    elif norm == 'pearson':
        if sp.sparse.issparse(D):
            D = D.A
            set_sparse = True
        else:
            set_sparse = False
        D = pearson_residuals(D, **args)
        if set_sparse:
            D = sp.sparse.csr_matrix(D)

    else:
        D.data[:] = (D.data / div)

    adata.X = D

def pearson_residuals(counts, theta=100):
    """
    Originally from https://github.com/berenslab/umi-normalization
    Computes analytical residuals for NB model with a fixed theta, clipping outlier residuals to sqrt(N)

    Parameters
    ----------
    counts : array-like
        Array of counts
    theta : float
        NB dispersion parameter

    Returns
    -------
    residuals : array-like
        Array of Pearson residuals
    """
    counts_sum0 = np.sum(counts, axis=0, keepdims=True)
    counts_sum1 = np.sum(counts, axis=1, keepdims=True)
    counts_sum  = np.sum(counts)

    #get residuals
    mu = counts_sum1 @ counts_sum0 / counts_sum
    z = (counts - mu) / np.sqrt(mu + mu**2/theta)

    #clip to sqrt(n)
    n = counts.shape[0]
    z[z >  np.sqrt(n)] =  np.sqrt(n)
    z[z < -np.sqrt(n)] = -np.sqrt(n)

    return z
