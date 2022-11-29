# Path: src/scripts/analysis.py
# Define functions used for analysis
# Author: Yuan Xue (xuesoso@gmail.com)
# Github: https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
# Updated on 28 Nov 2022
# License: MIT

#### Import variables and libraries ####
from ..variables.granuloma_lib import *
from .external import diffexp_10x
from .preprocess import *
from .utility import *
import collections
import h5py
from matplotlib.colors import to_hex
from glob import glob
import scvelo as scv


def run_scanpy(filename, adata=None, res=0.3, force=False, norm=None,
        log1p=True, include_genes=[]):
    """
    Run scanpy pipeline

    Parameters
    ----------
    filename : str
        Filename to save the output
    adata : anndata.AnnData
        AnnData object
    res : float
        Resolution for clustering
    force : bool
        Force to re-run the pipeline
    norm : str
        Normalization method
    log1p : bool
        Log1p transform the data
    include_genes : list
        List of genes to consider. If empty, run default variable gene filtering.

    Returns
    -------
    ad : anndata.AnnData
        Normalized AnnData object with PCA and UMAP objects computed.
    """
    adata_filename = filename + '.h5ad'
    if os.path.isfile(adata_filename) == False or force is True:
        print('Running Scanpy dimensionality reduction...')
        if len(include_genes) > 0:
            ad = adata.copy()[:, include_genes]
        else:
            ad = adata.copy()
        normalize_per_cell(ad, 'cell_median', max_fraction_genes=0.05)
        if norm is not None:
            normalize(ad, norm=norm)
        else:
            normalize(ad, norm='log1p')
        if norm == 'ftt':
            normalize(ad, norm='log1p')
        sc.pp.highly_variable_genes(ad, n_bins=40)
        ad_dr = sc.pp.filter_genes_dispersion(ad, copy=True)
        sc.pp.pca(ad_dr)
        sc.pp.neighbors(ad_dr, n_neighbors=50)
        sc.tl.umap(ad_dr)
        sc.tl.leiden(ad_dr, resolution=1)
        transfer_ad_key(ad.obs, ad_dr.obs, replace=True)
        transfer_ad_key(ad.uns, ad_dr.uns, replace=True)
        transfer_ad_key(ad.obsm, ad_dr.obsm, replace=True)
        transfer_ad_key(ad.obsp, ad_dr.obsp, replace=True)
        ad.uns['PCs'] = ad_dr.varm['PCs']
        ad.uns['neighbors']['connectivities'] = ad.obsp['connectivities']
        if 'dendrogram_leiden' in list(ad.uns.keys()):
            del(ad.uns['dendrogram_leiden'])
        sc.tl.rank_genes_groups(ad, 'leiden', method='wilcoxon')
        for key in ad.obs.columns:
            num = (len(sorted(set(ad.obs[key].values))))
            d = (ad.obs[key].values)
            if isinstance(d, pd.core.arrays.categorical.Categorical) and num == 1:
                ad.obs[key] = d.astype(str)
        save_h5ad(ad, adata_filename)
        del(ad_dr)
    else:
        print('file already exists. loading')
        ad = sc.readh5ad(adata_filename)
    return ad

def run_sam(filename, adata=None, force=False, res=0.3, **args):
    """
    A wrapper function to process scRNA-seq data with sam-algorithm (SAM).
    Refer to https://github.com/atarashansky/self-assembling-manifold for
    documentation on sam-algorithm.

    Writes SAM object, raw AnnData object, and processed AnnData object to the
    specified filename.

    Parameters
    ----------
    filename : str
        Filename to save the output
    adata : anndata.AnnData
        AnnData object
    force : bool
        Force to re-run the pipeline
    res : float
        Resolution for clustering
    **args : dict
        Additional arguments to pass to SAM.preprocess_data

    Returns
    -------
    out : anndata.AnnData
        AnnData object processed with SAM.
    """

    ## Running SAM
    if os.path.isfile(filename+'.p') is False or force is True:
        sam_wt = SAM(adata)
        print('file does not exist, loading...')
        sam_wt.preprocess_data(**args)
        sam_wt.run()
        sam_wt.leiden_clustering(res=res)
        m = sam_wt.identify_marker_genes_rf()
        sam_wt.save(filename)
        sam_wt.adata_raw.write_h5ad(filename+'_adata_raw.h5ad')
        sam_wt.adata.write_h5ad(filename+'_adata.h5ad')
        out =  sam_wt.adata
    else:
        print('file exists, loading...')
        out = anndata.readh5ad(filename+'_adata.h5ad')
    return out

def prepare_macrophage_marker(gpd, output_dir=infd+'../external/'):
    """
    Load canonical macrophage M1, M2, and lineage marker genes for visualization

    Parameters
    ----------
    gpd : dict
        A dictionary of gene name to gene ID dictionary
    output_dir : str
        Output directory

    Returns
    -------
    macrophage_marker : dict
        A dictionary of canonical macrophage marker genes
    """
    from matplotlib.colors import rgb2hex
    lineage_marker = np.intersect1d(list(gpd.keys()), [x.strip() for x in
                       pd.read_csv(output_dir+'lineage_marker.tsv',
                                   sep='   ')['Mouse'].values.flatten()])
    m1_marker = np.intersect1d(list(gpd.keys()), [x.strip() for x in
                                                  pd.read_csv(output_dir+'M1_marker.tsv',
                                                              sep='   ')['Mouse'].values.flatten()])
    m1_m2_marker = np.intersect1d(list(gpd.keys()), [x.strip()
                     for x in pd.read_csv(output_dir+'M1_M2_marker.tsv',
                                          sep='   ')['Mouse'].values.flatten()])
    m2_marker = np.intersect1d(list(gpd.keys()),
                               [x.strip() for x in pd.read_csv(output_dir+'M2_marker.tsv',
                                                               sep='   ')['Mouse'].values.flatten()])
    myeloid_marker = np.intersect1d(list(gpd.keys()),
                               [x.strip() for x in pd.read_csv(output_dir+'Myeloid_marker.tsv',
                                                               sep='   ')['Mouse'].values.flatten()])
    overlap_ind = np.in1d(list(gpd.keys()), np.concatenate([lineage_marker,
                                                            m1_marker, m1_m2_marker, m2_marker]))
    markers = np.array(list(gpd.keys()))[overlap_ind]
    df_markers = pd.DataFrame(markers, index=markers, columns=['marker'])
    df_markers['group'] = 'None'
    df_markers.loc[np.in1d(markers, m1_marker), 'group'] = 'M1'
    df_markers.loc[np.in1d(markers, m2_marker), 'group'] = 'M2'
    df_markers.loc[np.in1d(markers, np.setdiff1d(np.setdiff1d(m1_m2_marker,
                                              m2_marker), m1_marker)),
                   'group'] = 'M1|M2'
    df_markers.loc[np.in1d(markers, lineage_marker), 'group'] = 'lineage'
    df_markers = df_markers.sort_values('group')
    ucl = np.unique(df_markers.group)
    length =  len(ucl)
    return df_markers

def make_cmap(ad, groupby: str, cmap='tab20'):
    """
    Make a color map for a given groupby column in an AnnData object

    Parameters
    ----------
    ad : anndata.AnnData
        AnnData object
    groupby : str
        Column name in ad.obs to use for color map
    cmap : str
        Color map to use
    """
    assert groupby in ad.obs.columns, '{:} is not an available column'.format(ucl)
    ucl = sorted(set(ad.obs[groupby]))
    if type(cmap) is str:
        ad.uns[groupby+'_colors'] = np.array([to_hex(plt.get_cmap(cmap)(
            i/len(ucl))) for i, x in enumerate(ucl)])
    else:
        raise ValueError('Not implemented')

def run_velocyto(adh, fn, genes=None, force=False):
    """
    Run velocyto on the input anndata object.
    "adh" must contain neighbors and connectivities matrix

    Parameters
    ----------
    adh : anndata object
        anndata object with neighbors and connectivities matrix
    fn : str
        filename to save the output
    genes : list
        list of genes to use for velocity estimation
    force : bool
        force re-run of velocyto

    Returns
    -------
    adv : anndata object
        anndata object with RNA velocity information.
    """
    if os.path.isfile(fn) and force is False:
        print("Loading {:}".format(fn))
        adv = sc.readh5ad(fn)
    else:
        print("Running Velocyto analysis")
        adv = sc.readh5ad(infd+'201221_10X_velocyto_all.h5ad')
        adv.obs_names = [x.replace(':', '_').split('x-')[0]+'-1' for x in adv.obs_names]
        adv = adv[adh.obs_names]
        if genes is None:
            g = adh.var_names[adh.var['weights'] > 0]
        else:
            g = genes
        adv = adv[:, g]
        scv.pp.filter_and_normalize(adv)
        scv.pp.moments(adv)
        scv.tl.recover_dynamics(adv)
        scv.tl.velocity(adv, mode='dynamical')
        scv.tl.velocity_graph(adv)
        scv.tl.velocity_embedding(adv)
        transfer_ad_key(adv.obs, adh.obs)
        transfer_ad_key(adv.obsm, adh.obsm)
        adv.write_h5ad(fn)
    return adv

def run_velocyto_myeloid(adh, fn, genes=None, force=False, groupby='sub_leiden'):
    """
    Run velocyto on the myeloid cell subset.
    "adh" must contain neighbors and connectivities matrix

    Parameters
    ----------
    adh : anndata object
        anndata object with neighbors and connectivities matrix
    fn : str
        filename to save the output
    genes : list
        list of genes to use for velocity estimation
    force : bool
        force re-run of velocyto

    Returns
    -------
    adv : anndata object
        anndata object with RNA velocity information.
    """
    if os.path.isfile(fn) and force is False:
        print("Loading {:}".format(fn))
        adv = sc.readh5ad(fn)
    else:
        print("Running Velocyto analysis")
        adv = sc.readh5ad(infd+'201221_10X_velocyto_all.h5ad')
        adv.obs_names = [x.replace(':', '_').split('x-')[0]+'-1' for x in adv.obs_names]
        adv = adv[adh.obs_names]
        if genes is None:
            g = adh.var_names[adh.var['final_weights'] > 0]
        else:
            g = genes
        adv = adv[:, g]
        scv.pp.filter_and_normalize(adv)
        scv.pp.moments(adv)
        scv.tl.recover_dynamics(adv)
        scv.tl.velocity(adv, mode='dynamical')
        scv.tl.velocity_graph(adv)
        scv.tl.velocity_embedding(adv)
        transfer_ad_key(adv.obs, adh.obs)
        transfer_ad_key(adv.obsm, adh.obsm)
        sc.pl.umap(adh, color=groupby, show=False)
        adv.uns['{:}_colors'.format(groupby)] = adh.uns['{:}_colors'.format(groupby)].copy()
        scv.tl.paga(adv, groups=groupby, vkey='velocity')
        adv.write_h5ad(fn)
    return adv

def calculate_unit_norm(adh):
    """
    Calculate the unit norm for each cell

    Parameters
    ----------
    adh : AnnData
        AnnData object containing the data

    Returns
    -------
    unit_norm : np.array
        numpy array containing the normalized expression matrix
    """
    if sp.sparse.issparse(adh.X):
        unit_norm = (adh.X.A - np.mean(adh.X, axis=0).A) / np.std(adh.X.A, axis=0)
    else:
        unit_norm = (adh.X - np.mean(adh.X, axis=0)) / np.std(adh.X, axis=0)
    unit_norm = np.nan_to_num(unit_norm)
    return unit_norm

def calculate_correlation(unit_norm, target_vector):
    """
    Calculate the correlation between each cell and a target vector of cell
    attribute (e.g. latent biological pseudotime)

    Parameters
    ----------
    unit_norm : np.array
        numpy array containing the normalized expression matrix
    target_vector : np.array
        numpy array containing the target vector

    Returns
    -------
    df : pd.DataFrame
        pandas dataframe containing the correlation between each cell and the
        target vector
    """
    c1, c2, c3, c4, c5 = [], [], [], [], []
    r2 = []
    for x in unit_norm.T:
        slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x,
                target_vector)
        c1.append(slope)
        c2.append(intercept)
        c3.append(r_value)
        c4.append(p_value)
        c5.append(std_err)
        r2.append(r_value**2)
    df = pd.DataFrame(c1, index=range(len(unit_norm.T)),
            columns=['slope'])
    df['intercept'] = c2
    df['r_value'] = c3
    df['p'] = c4
    df['std_err'] = std_err
    df['r2'] = r2
    return df

def panglaoDB_score(adh, groupby='leiden', processed_panglao=None,
                    added_key='panglaoDB_cellType'):
    """
    Calculate the PanglaoDB cell type score for each cell. To do this, we compute
    a cell-type specific score for each cluster using a vector of cell-type
    specific genes and their associated weight. The weight vector down-weighs
    genes that are shared in expression across many cell types (non-uniqueness)
    and is provided by panglaoDB database (available in
    "data/external/PanglaoDB_markers_24_Jan_2020.tsv.gz").

    The full description of the method is available in the PanglaoDB paper.
    "Franzén et al., PanglaoDB: a web server for exploration of mouse and human
    single-cell RNA sequencing data, Database, Volume 2019, 2019, baz046,
    https://doi.org/10.1093/database/baz046"


    Parameters
    ----------
    adh : AnnData
        AnnData object containing the data
    groupby : str
        Key in adh.obs to group the samples by (e.g. clustering labels)
    processed_panglao : str
        Path to the processed panglaoDB marker gene table with weight vector
        (default: None)
    added_key : str
        Key to add the cell type score to in adh.obs (default: 'panglaoDB_cellType')

    Returns
    -------
    adh : AnnData
        AnnData object containing the data with cell type score
    """
    if processed_panglao is None:
        processed_panglao = infd+'../external/processed_PanglaoDB_markers_24_Jan_2020.tsv.gz'
    lineage_marker = load_panglaoDB(processed_panglao=processed_panglao)
    genes = lineage_marker['GeneSymbol'].values
    groups = adh.obs[groupby].values
    unit_norm = calculate_unit_norm(adh)
    unique_lineage = sorted(set(lineage_marker['cell type']))
    unique_groups = sorted(set(groups))
    df_lineage = pd.DataFrame(0, index=unique_lineage, columns=unique_groups)
    max_frequency = max(lineage_marker['frequency'].values)
    min_frequency = min(lineage_marker['frequency'].values)
    __range__ = np.arange(len(adh.var_names))
    for c in unique_groups:
        subset = (groups == c)
        for i in unique_lineage:
            row_subset = (lineage_marker['cell type'].values == i)
            subgenes = sorted(set(genes[row_subset]))
            mask = np.squeeze([__range__[adh.var_names == x] for x in subgenes])
            median_exp = np.mean(unit_norm[:, mask][subset], axis=0)
            weight = 1 + np.sqrt((np.max(max_frequency - lineage_marker.loc[
                row_subset, 'frequency'].values))/(max_frequency - min_frequency))
            weighted_score = np.sum(median_exp * weight) / len(subgenes)**(1/3)
            df_lineage.loc[i, c] = weighted_score
    key, values = df_lineage.idxmax().reset_index().values.T
    winner = pd.DataFrame([key, values], index=['cluster', 'cell_type'], columns=key).T
    specificity = []
    for i in winner.cluster:
        score = (df_lineage[winner['cluster'][i]].sort_values())
        diff = score[-1] - score[-2]
        score_range = max(score) - min(score)
        specificity.append(diff/score_range)
    df = pd.DataFrame(specificity, index=winner.cluster, columns=['specificity'])
    adh.obs[added_key] = [winner.loc[x, 'cell_type'] for x in adh.obs[groupby]]
    df_lineage.columns = df_lineage.columns.values.astype(str)
    adh.uns[added_key] = df_lineage

#### Differential expression functions

def multiple_test_correction(sorted_pvals):
    """
    calculate false discovery rate (FDR) by performing Benjamin-Hochberg
    correction.

    Parameters
    ----------
    sorted_pvals : list
                   Sorted by ascending order p-values.

    Returns
    -------
    adjusted : np.array
            Sorted adjusted p-values (FDR).
    """

    ## Benjamin-Hochberg correction
    n = np.shape(sorted_pvals)[0]
    m = np.arange(1, n+1)
    adjusted = sorted_pvals*n/m
    adjusted[adjusted > 1] = 1
    return adjusted

def enrich_genes(exp, C, groupby='specific', method='wilcoxon'):
    """
    Identify differential gene expression between two groups.

    Parameters
    ----------
    exp : pandas.DataFrame
            Expression matrix with genes as index and cells as columns.
    C : pandas.Series
            Cell labels.
    groupby : str
            Grouping variable. Default is 'specific'.
    method : str
            Choice of statistic test. Default is 'wilcoxon'. Other options
            are 'ttest' and 'mannwhitneyu'.

    Returns
    -------
    out : pandas.DataFrame
            Table with differential expression results.
    """
    if isinstance(exp, anndata.AnnData):
        exp, C, __ = from_adata_to_df(exp, obs=C)
    if groupby=='specific':
        UC = np.unique(C)
        avg = average_by_clusters(exp, C, return_clusters=True)
        UC_array = np.repeat(UC, avg.shape[1]).reshape(avg.shape[0], -1)
        id_index = avg.index.values
        length = len(UC)
        max_exp_cluster = np.array([UC[x] for x in np.argmax(avg.values, axis=0)])
        max2_exp_cluster = np.array([UC[x] for x in
                                     np.argsort(avg.values, axis=0)[-2]])
        E = exp.values
        outvals = []
        for eid, g in enumerate(exp.columns.values):
            sub1 = max_exp_cluster[eid]
            sub2 = max2_exp_cluster[eid]
            v1 = E[C == sub1, eid]
            v2 = E[C == sub2, eid]
            if method == 'wilcoxon':
                _stats, p_val = sp.stats.ranksums(v1, v2)
            elif method == 'ttest':
                _stats, p_val = sp.stats.ttest_ind(v1, v2)
            elif method == 'mannwhitneyu':
                try:
                    _stats, p_val = sp.stats.mannwhitneyu(v1, v2)
                except:
                    _stats = 0
                    p_val = 1
            else:
                raise ValueError('Invalid method')
            outvals.append([v1.mean()/v2.mean(), v1.mean(), v2.mean(), p_val, sub1, sub2])
        out = pd.DataFrame(outvals, index=exp.columns, columns=['fold', 'max1_mean', 'max2_mean', 'p',
                                                                'max1', 'max2'])
        out = out.fillna(0)
    else:
        UC = np.unique(C)
        avg = average_by_clusters(exp, C, return_clusters=True)
        UC_array = np.repeat(UC, avg.shape[1]).reshape(avg.shape[0], -1)
        id_index = avg.index.values
        length = len(UC)
        max_exp_cluster = np.array([UC[x] for x in np.argmax(avg.values, axis=0)])
        E = exp.values
        outvals = []
        for eid, g in enumerate(exp.columns.values):
            sub1 = max_exp_cluster[eid]
            v1 = E[C == sub1, eid]
            v2 = E[C != sub1, eid]
            _stats, p_val = sp.stats.ranksums(v1, v2)
            if method == 'wilcoxon':
                _stats, p_val = sp.stats.ranksums(v1, v2)
            elif method == 'ttest':
                _stats, p_val = sp.stats.ttest_ind(v1, v2)
            elif method == 'mannwhitneyu':
                try:
                    _stats, p_val = sp.stats.mannwhitneyu(v1, v2)
                except:
                    _stats = 0
                    p_val = 1
            else:
                raise ValueError('Invalid method')
            outvals.append([v1.mean()/v2.mean(), p_val, sub1])
        out = pd.DataFrame(outvals, index=exp.columns, columns=['fold', 'p', 'max1'])
        out = out.fillna(0)
    out['log2-fold'] = np.log2(out['fold'].values)
    noFoldChange = (out['fold'].values == 1)
    out.loc[noFoldChange, 'log2-fold'] = 0
    out['abs-log2-fold'] = np.abs(out['log2-fold'].values)
    out = out.sort_values(['p', 'abs-log2-fold'], ascending=[True, False])
    out['adjusted-p'] = multiple_test_correction(out['p'].values)
    out.loc[out['adjusted-p'] > 1] = 1
    return out

def export_marker(ad, directory, groupby='leiden'):
    """
    A wrapper function to export marker genes for each cluster to
    "{directory}/{cluster}/".

    Parameters
    ----------
    ad : :class:`~anndata.AnnData`
        Annotated data matrix.
    directory : str
                Directory to save the marker genes.
    groupby : str
                Key in `ad.obs` to groupby.
    """
    os.makedirs(directory, exist_ok=True)
    if sp.sparse.issparse(ad.X):
        D = pd.DataFrame(ad.X.A, index=ad.obs_names, columns=ad.var_names)
    else:
        D = pd.DataFrame(ad.X, index=ad.obs_names, columns=ad.var_names)
    for f in glob(directory+'/*/*.tsv.gz'):
        os.remove(f)
    os.makedirs(directory+'specific_marker/', exist_ok=True)
    marker = enrich_genes(D, ad.obs[groupby].values,
            groupby='specific')
    for i in sorted(set(ad.obs[groupby].values)):
        subset = marker[marker['max1'] == i]
        print(i)
        subset.to_csv(directory+'specific_marker/cluster_'+str(i).replace('/', '_')+'_marker.tsv.gz', sep='\t', compression='gzip')
    os.makedirs(directory+'negative_binomial/', exist_ok=True)
    dtest = diffexp_10x.run(ad, groupby, layer='raw_X', directory=directory+'negative_binomial/')
    print('Exported marker genes!')

def cmh_test(df, K, verbose=False):
    """
    Performs and retrieves Cochran-Mantel-Haenszel test statistic
    input Table should be arranged in pairs.
    for instance, we have a table

        A1 A2 A3 A4 B1 B2 B3 B4
    X   .. .. .. .. .. .. .. ..
    Y   .. .. .. .. .. .. .. ..
    Z   .. .. .. .. .. .. .. ..


    where A1, A2, A3, A4 are the cell counts from mice treated with condition A
    and B1, B2, B3, B4 are the cell counts from mice treated with condition B.
    If we set K as 4, then we will pair A1 with B1, A2 with B2, and so on...
    We then perform CMH chi-squared test to assess conditional independence.

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas DataFrame containing the cell counts.
    K : int
        Number of pairs to be tested.
    verbose : bool
        Whether to print the test results.

    Returns
    -------
    cmh_res : dict.
        A dictionary containing the Chi-squared test model and the test statistic.
    cmh_table : pandas.DataFrame
        A pandas DataFrame containing the test results.
    """

    #### Convert to 2 x 2 x K contignecy table
    size = df.sum(0)
    cmh_res = {}
    for x in df.index:
        tables = []
        for i in range(K):
            x1, x2 = df.loc[x].values[i], df.loc[x].values[i+K]
            n1, n2 = size[i], size[i+K]
            table = [[x1, n1-x1], [x2, n2-x2]]
            tables.append(table)
        st = sm.stats.StratifiedTable(tables)
        cmhstat, pval = str(st.summary()).split('\nTest of OR=1')[-1].split('\n')[0].strip().split('   ')
        st.stat = float(cmhstat)
        st.pval = float(pval)
        cmh_res[x] = st

    ####  Multiple hypothesis adjustment
    cmh_table = pd.DataFrame(0, index=df.index, columns=['cmh.stat', 'pval', 'fdr',
        'oddsratio_pooled', 'riskratio_pooled', 'logodds_pooled', 'logodds_pooled_se',
        'oddsratio_pooled_confint_min', 'oddsratio_pooled_confint_max'])

    pvals = []
    for x in cmh_res:
        pvals.append(cmh_res[x].pval)
    pvals = np.array(pvals)
    order = np.argsort(pvals)
    fdr = multiple_test_correction(pvals[order])

    for i, x in enumerate(np.array(list(cmh_res.keys()))[order]):
        mod = cmh_res[x]
        mod.fdr = fdr[i]
        cmh_table.loc[x, 'cmh.stat'] = mod.stat
        cmh_table.loc[x, 'pval'] = mod.pval
        cmh_table.loc[x, 'fdr'] = mod.fdr
        cmh_table.loc[x, 'ratio'] = mod.oddsratio_pooled
        cmh_table.loc[x, 'log2_ratio'] = np.log2(mod.oddsratio_pooled)
        cmh_table.loc[x, 'oddsratio_pooled'] = mod.oddsratio_pooled
        cmh_table.loc[x, 'riskratio_pooled'] = mod.riskratio_pooled
        cmh_table.loc[x, 'logodds_pooled'] = mod.logodds_pooled
        cmh_table.loc[x, 'logodds_pooled_se'] = mod.logodds_pooled_se
        oddsratio_pooled_confint = mod.oddsratio_pooled_confint()
        cmh_table.loc[x, 'oddsratio_pooled_confint_min'] = oddsratio_pooled_confint[0]
        cmh_table.loc[x, 'oddsratio_pooled_confint_max'] = oddsratio_pooled_confint[1]

    if verbose:
        print(cmh_table)
    return cmh_res, cmh_table

def run_GSOA(fns, outDir, groups, force=False, verbose=False):
    """
    Load differentially expressed marker genes and compute gene set
    over representation analysis (GSOA).

    Parameters
    ----------
    fns : list
        List of files containing negative binomial differential genes.
    outDir : str
        Output directory.
    groups : list
        List of group factors (i.e. cluster names).
    force : bool
        Force re-run of GSOA.
    verbose : bool
        Whether to print the output directory path

    Returns
    -------
    load_GSOA(outDir) : dict
        A dictionary containing the GSOA results.
    """
    if verbose:
        print('designated output directory: {:}'.format(outDir))
    import shutil

    hm = read_msigdb(infd+'../external/MSigDB/h.all.v7.4.symbols.gmt')
    tft = read_msigdb(infd+'../external/MSigDB/c3.tft.v7.4.symbols.gmt')
    imgs = read_msigdb(infd+'../external/MSigDB/c7.immunesigdb.v7.4.symbols.gmt')
    kegg_metabolism = read_msigdb(infd+'../external/KEGG/KEGG_mmu01100_pathways.txt')
    biocarta = read_msigdb(infd+'../external/MSigDB/c2.cp.biocarta.v7.4.symbols.gmt')
    pid = read_msigdb(infd+'../external/MSigDB/c2.cp.pid.v7.4.symbols.gmt')
    wikipathways = read_msigdb(infd+'../external/MSigDB/c2.cp.wikipathways.v7.4.symbols.gmt')
    reactome = read_msigdb(infd+'../external/MSigDB/c2.cp.reactome.v7.4.symbols.gmt')
    os.makedirs(outDir, exist_ok=True)

    l2thres, pthres, min_expressed = 0.25, 0.05, 0.25

    marker, marker_down = {}, {}
    for fn in fns:
        cl = fn.split('/')[-1].split('_')[0]
        df = pd.read_csv(fn, sep='\t', index_col=0)
        marker[cl] = df[(df['log2-fold']>l2thres)*(df['adjusted-p']<pthres)*(df['expressed_a']>min_expressed)].sort_values('log2-fold', ascending=False)
        marker_down[cl] = df[(df['log2-fold']<-l2thres)*(df['adjusted-p']<pthres)*(df['expressed_b']>min_expressed)].sort_values('log2-fold', ascending=True)

    if os.path.isfile(outDir+'{:}_HallmarkUp.tsv.gz'.format(groups[0])) is False or force:
        if os.path.isdir(outDir) and outDir.split('/')[-1] == 'MSigDB':
            shutil.rmtree(outDir)
        os.makedirs(outDir, exist_ok=True)
        hmres = hypergeometric_test(marker, groups, 'HallmarkUp', hm, outDir)
        tftres = hypergeometric_test(marker, groups, 'TargetTranscriptionUp', tft, outDir)
        imgres = hypergeometric_test(marker, groups, 'ImmuneSignatureUp', imgs, outDir)
        metabolism_res = hypergeometric_test(marker, groups, 'KEGG_MetabolismUp', kegg_metabolism, outDir)
        biocartares = hypergeometric_test(marker, groups, 'BioCartaUp', biocarta, outDir)
        pidres = hypergeometric_test(marker, groups, 'PID_Up', pid, outDir)
        wikipathwaysres = hypergeometric_test(marker, groups, 'WikipathwaysUp', wikipathways, outDir)
        reactomeres = hypergeometric_test(marker, groups, 'ReactomeUp', reactome, outDir)

        hmres_down = hypergeometric_test(marker_down, groups, 'HallmarkDown', hm, outDir)
        tftres_down = hypergeometric_test(marker_down, groups, 'TargetTranscriptionDown', tft, outDir)
        imgres_down = hypergeometric_test(marker_down, groups, 'ImmuneSignatureDown', imgs, outDir)
        metabolism_res_down = hypergeometric_test(marker_down, groups, 'KEGG_MetabolismDown', kegg_metabolism, outDir)
        biocartares_down = hypergeometric_test(marker_down, groups, 'BioCartaDown', biocarta, outDir)
        pidres_down = hypergeometric_test(marker, groups, 'PID_Down', pid, outDir)
        wikipathwaysres_down = hypergeometric_test(marker_down, groups, 'WikipathwaysDown', wikipathways, outDir)
        reactomeres_down = hypergeometric_test(marker_down, groups, 'ReactomeDown', reactome, outDir)

    return load_GSOA(outDir)

def _hypergeometric_test(query, gene_sets, N_genes=31053):
    """
    Perform hypergeometric test to assess over-representation of query genes
    in a reference gene set.

    Parameters
    ----------
    query : the list of genes to test for
    gene_sets : the dictionary of gene_sets with keys being the pathway name
    N_genes : the total gene set space (31053 is the number for Mus_musculus 10)

    Returns
    -------
    result : a dataframe with the results of the test
    """
    result = {}
    for key in gene_sets:
        set1 = query
        set2 = gene_sets[key]
        int_set = np.intersect1d(set1, set2)

        N = N_genes
        n = len(set1)
        K = len(set2)
        k = len(int_set)

        pval = sp.stats.hypergeom(N, K, n).sf(k-1)
        result[key] = [key, k, n, K, ','.join(int_set), pval]
    result = pd.DataFrame(result, index=['pathway', 'intersect_size', 'query_size', 'set_size', 'intersect_gene', 'pval']).T.sort_values('pval')
    result['padj'] = result['pval']*np.arange(1, len(result)+1)
    result.loc[result['padj'] > 1, 'padj'] = 1
    return result

def hypergeometric_test(marker, keys, name, genesets, outDir):
    """
    Wrapper function to perform hypergeometric test and save the results of the
    pathway test to a file.

    Parameters
    ----------
    marker : the dictionary of marker genes with keys being the cluster name
    keys : the list of cluster names
    name : the name of the test
    genesets : the dictionary of gene_sets with keys being the pathway name
    outDir : the directory to save the results

    Returns
    -------
    pathwayres : a dataframe with the test results
    """
    pathwayres = {}
    for key in keys:
        filename = outDir + '{:}_{:}.tsv.gz'.format(key, name)
        res = _hypergeometric_test(marker[key].index.values, genesets)
        pathwayres[key] = res
        res.to_csv(filename, sep='\t', compression='gzip')
    return pathwayres
