# Path: src/visualization/plot.py
# A library of functions for visualizing scRNA-seq data analysis.
# Author: Yuan Xue (xuesoso@gmail.com)
# Github: https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
# Updated on 28 Nov 2022
# License: MIT

#### Import libraries ####
import warnings
warnings.filterwarnings("ignore")
import sys, os, subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib import gridspec
from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import ListedColormap, Normalize
import scipy as sp
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
from DensityPlot import density2d
import pandas.io.formats.style
import anndata
from adjustText import adjust_text
import seaborn as sns


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

class MathTextSciFormatter(mticker.Formatter):
    '''
    This formatter can be fed to set ticklabels in scientific notation without
    the annoying "1e" notation (why would anyone do that?).
    Instead, it formats ticks in proper notation with "10^x".

    fmt: the decimal point to keep
    Usage = ax.yaxis.set_major_formatter(MathTextSciFormatter("%1.2e"))
    '''

    def __init__(self, fmt="%1.2e"):
        self.fmt = fmt
    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        tup = s.split('e')
        significand = tup[0].rstrip(decimal_point)
        sign = tup[1][0].replace(positive_sign, '')
        exponent = tup[1][1:].lstrip('0')
        if exponent:
            exponent = '10^{%s%s}' % (sign, exponent)
        if significand and exponent:
            s =  r'%s{\times}%s' % (significand, exponent)
        else:
            s =  r'%s%s' % (significand, exponent)
        return "${}$".format(s)

class FormatScalarFormatter(mticker.ScalarFormatter):
    def __init__(self, fformat="%1.1f", offset=True, mathText=True):
        self.fformat = fformat
        mticker.ScalarFormatter.__init__(self, useOffset=offset,
                                            useMathText=mathText)

    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % mticker._mathdefault(self.format)

class scatter():
    def __init__(self,
                 X,
                 c=None,
                 filename='',
                 figsize=(3, 3),
                 cmap=None,
                 s=5,
                 title=True,
                 legend_title=None,
                 facecolor='white',
                 xlabel=None,
                 ylabel=None,
                 show_cbar=True,
                 alpha=1,
                 fontsize=13,
                 hspace=0.15,
                 readable_titles=False,
                 wspace=0.05,
                 insertion_spacer=0,
                 clusters=None,
                 rotation=0,
                 ncols=5,
                 sort_by_str=None,
                 sina=False,
                 run=True,
                 cbar_ticks=True,
                 density=False,
                 edgecolor=None,
                 linewidths=None,
                 density_cmap='Reds',
                 bw='scott',
                 density_alpha=0.6,
                 dot_order='ordered',
                 inverse_dot_order=False,
                 order_cbar_labels=[],
                 adj=None,
                 showxticks=True,
                 showyticks=True,
                 ylim=None,
                 xlim=None,
                 legend_loc=None,
                 vmin=None,
                 vmax=None,
                 stratify=False,
                 grid=False,
                 labelsize=12,
                 markerscale=None,
                 bins=300,
                 axis_border='on',
                 force_discrete_cbar=False,
                 sharex=True,
                 sharey=True,
                 fn='',
                 layer=None,
                 symbol=None,
                 components=[0, 1],
                 dpi=300,
                 **args):
        self.label_dtype = []
        if isinstance(X, anndata.AnnData):
            self._is_anndata = True
            if layer is not None:
                self.X = X.obsm['X_'+layer][:, components]
            else:
                if 'X_umap' in X.obsm.keys():
                    self.X = X.obsm['X_umap'][:, components]
                    layer = 'UMAP'
                elif 'X_pca' in X.obsm.keys():
                    self.X = X.obsm['X_pca'][:, components]
                    layer = 'PCA'
                else:
                    raise ValueError('No valid layer is provided')
            if isinstance(c, (str, list, pd.core.series.Series, np.ndarray, pd.core.indexes.base.Index)):
                if title:
                    if symbol is not None:
                        title = X.var.loc[c, symbol]
                    elif stratify:
                        title = sorted(set(X.obs[c].values))
                    elif c is not None:
                        title = c
                    else:
                        title = ''
                if type(c) is str:
                    c = [c]
                _c = []
                for _plot_num, ci in enumerate(c):
                    assert ci in X.obs.columns or ci in X.var_names, '{:} is not a valid observation or feature key'.format(c)
                    if ci in X.obs.columns:
                        _cvals = X.obs[ci].values
                        _c.append(_cvals)
                        _c_dtype = X.obs[ci].values.dtype
                        self.label_dtype.append(_c_dtype)
                    else:
                        if sp.sparse.issparse(X.X):
                            _cvals = X[:, ci].X.A.flatten()
                        else:
                            _cvals = X[:, ci].X.flatten()
                        _c.append(_cvals)
                        _c_dtype = 'float'
                        self.label_dtype.append(_c_dtype)
                if len(c) <= 1:
                    c = np.array(_c[0])
                else:
                    c = np.array(_c)
            if type(adj) is str:
                adj = X.uns[adj]
            elif adj is True:
                if 'connectivities' in X.obsp.keys():
                    adj = X.obsp['connectivities']
                else:
                    adj = X.uns['neighbors']['connectivities']
            if type(clusters) is str:
                clusters = X.obs[clusters].values
        else:
            self._is_anndata = False
            self.X = X
        if c is None:
            self.c = np.zeros(np.shape(self.X)[0])
            show_cbar = False
        else:
            ndim_c = np.ndim(c)
            if adj is not None:
                assert np.shape(adj)[0] == np.shape(adj)[1],\
                        'Adjacency matrix should have syemmetric length'
                assert np.shape(adj)[0] == np.shape(self.X)[0],\
                        'Adjacency matrix should have the same length as\
                        number of samples'
                if ndim_c == 1:
                    self.c = np.array(adj.dot(c) / np.sum(adj, axis=1).reshape(
                        np.shape(c)[0], 1).T)[0]
                else:
                    self.c = average_by_nnm(c.T, adj).T
            else:
                if isinstance(c, pd.Categorical):
                    c = c.astype(str)
                self.c = np.array(c)
                if self._is_anndata is False:
                    for ci in self.c:
                        self.label_dtype.append(ci.dtype)
        self.clusters = clusters
        if self.clusters is not None:
            assert np.shape(clusters)[0] == np.shape(self.X)[0],\
                    'the length of clusters must equal the length of input\
                    samples'
            figsize = (figsize[0], figsize[1] + figsize[0] / 4)
            self.clusters = np.array(clusters)
        self.show_cbar = show_cbar
        if fn != '':
            self.filename = fn
        else:
            self.filename = filename
        self.legend_loc = legend_loc
        self.legend_title = legend_title
        if self.filename != '':
            self.save_plot = True
            dirname = os.path.dirname(os.path.abspath(self.filename))
            if os.path.isdir(dirname) == False:
                os.mkdir(dirname)
        else:
            self.save_plot = False
        self.column_size = int(ncols)
        self.figsize = figsize
        if stratify == True and cmap is None:
            self.cmap = {x:'C'+str(y) for y, x in enumerate(sorted(set(self.c)))}
        else:
            self.cmap = cmap
        if type(c) is not str and type(title) is bool:
            self.title = ''
        else:
            self.title = title
        self.facecolor = facecolor
        if xlabel or xlabel == '':
            self.xlabel = xlabel
        else:
            if layer:
                self.xlabel = layer.split('_')[-1].upper()+str(components[0]+1)
            else:
                self.xlabel = None
        if ylabel or xlabel == '':
            self.ylabel = ylabel
        else:
            if layer:
                self.ylabel = layer.upper().split('_')[-1]+str(components[1]+1)
            else:
                self.ylabel = None
        self.alpha = alpha
        self.fontsize = fontsize
        self.insertion_spacer = insertion_spacer
        self.grid = grid
        self.markerscale = markerscale
        self.labelsize = labelsize
        if self.insertion_spacer > 0:
            self.title = make_titles_readable(self.title,
                                    insertion_spacer=self.insertion_spacer)
        self.s = s
        self.hspace = hspace
        self.rotation = rotation
        self.sina = sina
        self.cbar_ticks = cbar_ticks
        self.density = density
        self.bw = bw
        self.density_alpha = density_alpha
        self.fig_transparent = False
        self.dot_order = dot_order
        self.inverse_dot_order = inverse_dot_order
        self.order_cbar_labels = order_cbar_labels
        self.edgecolor = edgecolor
        self.linewidths = linewidths
        self.showxticks = showxticks
        self.showyticks = showyticks
        self.bins = bins
        self.args = args
        self.ylim = ylim
        self.xlim = xlim
        self.vmin = vmin
        self.vmax = vmax
        self.stratify = stratify
        self.axis_border = axis_border
        self.force_discrete_cbar = force_discrete_cbar
        self.sharex = sharex
        self.sharey = sharey
        self.dpi = dpi
        if self.facecolor == 'w' or self.facecolor == 'white':
            self.fig_transparent = True
        if self.density == False and density_cmap == 'None':
            self.density_cmap = self.cmap
        else:
            self.density_cmap = density_cmap
        if sort_by_str is not None and self.clusters is not None:
            assert np.in1d(sort_by_str, np.unique(self.clusters)).sum() ==\
                    len(np.unique(self.clusters)), 'The provided sort_by_str\
                    must contain all unique cluster values'

            self.sort_by_str = np.array(sort_by_str)
        else:
            self.sort_by_str = sort_by_str
        if self.show_cbar == True:
            self.wspace = (wspace + 0.15)
        else:
            self.wspace = wspace
        self.readable_titles = readable_titles
        plt.rcParams["axes.edgecolor"] = "0.15"
        plt.rcParams["axes.linewidth"] = 1.25
        if run == True:
            self.run()

    def run(self):
        plot_ndim = self.c.ndim
        if plot_ndim == 1:
            if self.force_discrete_cbar == False:
                self.color_dtype = self.c.dtype
            else:
                self.color_dtype = 'str'
            if self.stratify == False:
                return (self.make_single_plot(X=self.X, c=self.c,
                                              save_plot=self.save_plot))
            else:
                return (self.make_multiple_plot(save_plot=self.save_plot))
        else:
            self.color_dtype = self.label_dtype
            return (self.make_multiple_plot(save_plot=self.save_plot))

    def make_single_plot(self, X, c, save_plot=False, legend_loc=None):
        fig, ax = plt.subplots(figsize=self.figsize)
        self.figure = fig
        self.ax = ax
        if self.clusters is not None:
            gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])
            ax_bottom = plt.subplot(gs[1])
            ax = plt.subplot(gs[0])
            ## We get the old position and shift the bottom plot a bit up to
            ## remove dead space and give more space for the titles
            curr_pos = ax_bottom.get_position()
            curr_pos_top = ax.get_position()
            average_height = (curr_pos.height + curr_pos_top.height) / 2
            ## We set a new position for the bottom plot so that it has more
            ## space between the current top and next top plot
            new_pos = [
                curr_pos.x0,
                (curr_pos.y0 + (self.hspace - 0.05) / 2 * curr_pos.height +
                 (self.hspace - 0.05) / 2 * curr_pos_top.height),
                curr_pos.width, curr_pos.height
            ]
            ax_bottom.set_position(new_pos)
            self.make_bottom_plot(ax_bottom=ax_bottom)
        if self.ylim:
            ax.set_ylim(self.ylim)
        if self.xlim:
            ax.set_xlim(self.xlim)
        scatter_plot(X=X,
                     ax=ax,
                     fig=fig,
                     c=c,
                     xlabel=self.xlabel,
                     ylabel=self.ylabel,
                     fontsize=self.fontsize,
                     title=self.title,
                     legend_title=self.legend_title,
                     number_of_plots=1,
                     facecolor=self.facecolor,
                     cmap=self.cmap,
                     s=self.s,
                     alpha=self.alpha,
                     edgecolor=self.edgecolor,
                     linewidths=self.linewidths,
                     show_cbar=self.show_cbar,
                     column_size=self.column_size,
                     cbar_ticks=self.cbar_ticks,
                     density=self.density,
                     density_cmap=self.density_cmap,
                     bw=self.bw,
                     density_alpha=self.density_alpha,
                     dot_order=self.dot_order,
                     inverse_dot_order=self.inverse_dot_order,
                     order_cbar_labels=self.order_cbar_labels,
                     showxticks=self.showxticks,
                     showyticks=self.showyticks,
                     legend_loc=self.legend_loc,
                     vmin=self.vmin,
                     vmax=self.vmax,
                     grid=self.grid,
                     markerscale=self.markerscale,
                     labelsize=self.labelsize,
                     bins=self.bins,
                     axis_border=self.axis_border,
                     color_dtype=self.color_dtype,
                     **self.args)
        if save_plot != False and self.show_cbar == True:
            fig.savefig(self.filename,
                        bbox_inches='tight',
                        dpi=self.dpi,
                        transparent=self.fig_transparent)
        elif save_plot != False and self.show_cbar == False:
            fig.savefig(self.filename,
                        bbox_inches='tight',
                        dpi=self.dpi,
                        transparent=self.fig_transparent)
        else:
            return (fig, ax)

    def make_multiple_plot(self, save_plot=False):
        if self.stratify == False:
            number_of_plots = np.shape(self.c)[0]
        else:
            assert self.color_dtype not in ['int32', 'int64', 'float'],\
                    'must provide categorical color types'
            number_of_plots = len(np.unique(self.c))
        print('Number of plots: %s' % number_of_plots)
        data_array = np.array(self.X)
        if type(self.title) == str:
            new_title = np.repeat(self.title, number_of_plots)
        else:
            new_title = self.title
        assert np.shape(new_title)[0] == number_of_plots,\
        'the number of provided titles does not match that of the\
                provided colors'
        if self.readable_titles:
            new_title = make_titles_readable(titles=new_title,
                                             insertion_spacer=\
                                             self.insertion_spacer)
        row_size = int(np.ceil(number_of_plots / self.column_size))
        spacer_columns = 0
        if self.show_cbar:
            new_figsize = ((self.figsize[0]) * (self.column_size +
                spacer_columns) * 1.2, self.figsize[1] * row_size)
        else:
            new_figsize = ((self.figsize[0]) * (self.column_size +
                spacer_columns), self.figsize[1] * row_size)
        if self.clusters is None:
            fig, axs = plt.subplots(ncols=self.column_size,
                                    nrows=row_size,
                                    sharex=self.sharex,
                                    sharey=self.sharey,
                                    figsize=(new_figsize))
            axs = axs.flatten()
            for ax in axs[number_of_plots:]:
                ax.set_visible(False)
            fig.subplots_adjust(wspace=self.wspace, hspace=self.hspace)
        else:
            widths = np.repeat(4, self.column_size)
            heights = np.tile([4, 1], row_size)
            self.hspace *= 2
            fig = plt.figure(figsize=(new_figsize))
            new_grids = gridspec.GridSpec(ncols=self.column_size,
                                          nrows=row_size * 2,
                                          wspace=self.wspace,
                                          width_ratios=widths,
                                          height_ratios=heights,
                                          hspace=self.hspace)
            axs = [fig.add_subplot(x) for x in new_grids]
        self.figure = fig
        self.ax = axs
        count = 0
        row_count_1 = 0
        row_count_2 = self.column_size
        user_supplied_label_dtype = (len(self.label_dtype) > 0)
        color_dtype_ndim = np.ndim(self.color_dtype)
        curr_color_dtype = self.color_dtype
        if self.stratify == True:
            to_color = sorted(set(self.c))
        else:
            to_color = self.c
        for c, title in zip(to_color, new_title):
            if self.stratify == True:
                sub = (self.c == c)
                X = data_array[sub]
                c = np.repeat(c, sub.sum())
            else:
                X = data_array
            if count < number_of_plots:
                if self.clusters is None:
                    ax = axs[count]
                    if self.ylim:
                        ax.set_ylim(self.ylim)
                    if self.xlim:
                        ax.set_xlim(self.xlim)
                    if color_dtype_ndim > 0:
                        if self.color_dtype[count] in ['float', 'int64', 'int32', 'float32', 'float64']:
                            c = c.astype(float)
                        else:
                            c = c.astype(str)
                        curr_color_dtype = self.color_dtype[count]
                    scatter_plot(X=X,
                                 ax=ax,
                                 fig=fig,
                                 c=c,
                                 xlabel=self.xlabel,
                                 ylabel=self.ylabel,
                                 fontsize=self.fontsize,
                                 title=title,
                                 legend_title=self.legend_title,
                                 number_of_plots=number_of_plots,
                                 facecolor=self.facecolor,
                                 cmap=self.cmap,
                                 s=self.s,
                                 alpha=self.alpha,
                                 edgecolor=self.edgecolor,
                                 linewidths=self.linewidths,
                                 show_cbar=self.show_cbar,
                                 column_size=self.column_size,
                                 cbar_ticks=self.cbar_ticks,
                                 density=self.density,
                                 density_cmap=self.density_cmap,
                                 bw=self.bw,
                                 density_alpha=self.density_alpha,
                                 dot_order=self.dot_order,
                                 inverse_dot_order=self.inverse_dot_order,
                                 order_cbar_labels=self.order_cbar_labels,
                                 vmin=self.vmin,
                                 vmax=self.vmax,
                                 grid=self.grid,
                                 markerscale=self.markerscale,
                                 labelsize=self.labelsize,
                                 bins=self.bins,
                                 axis_border=self.axis_border,
                                 showxticks=self.showxticks,
                                 showyticks=self.showyticks,
                                 legend_loc=self.legend_loc,
                                 color_dtype=curr_color_dtype)
                else:
                    ax0 = axs[row_count_1]
                    ax1 = axs[row_count_2]
                    ## We get the old position and shift the bottom plot a bit
                    ## up to remove dead
                    ## space and give more space for the titles
                    curr_pos = ax1.get_position()
                    curr_pos_top = ax0.get_position()
                    average_height = (curr_pos.height +
                                      curr_pos_top.height) / 2
                    ## We set a new position for the bottom plot so that it
                    ## has more space between the current top and next top plot
                    new_pos = [
                        curr_pos.x0,
                        (curr_pos.y0 +
                         (self.hspace - 0.05) / 2 * curr_pos.height +
                         (self.hspace - 0.05) / 2 * curr_pos_top.height),
                        curr_pos.width, curr_pos.height
                    ]
                    ax1.set_position(new_pos)
                    scatter_plot(X=X,
                                 ax=ax0,
                                 fig=fig,
                                 c=c,
                                 xlabel=self.xlabel,
                                 ylabel=self.ylabel,
                                 fontsize=self.fontsize,
                                 title=title,
                                 legend_title=self.legend_title,
                                 number_of_plots=number_of_plots,
                                 facecolor=self.facecolor,
                                 cmap=self.cmap,
                                 s=self.s,
                                 alpha=self.alpha,
                                 edgecolor=self.edgecolor,
                                 linewidths=self.linewidths,
                                 show_cbar=self.show_cbar,
                                 column_size=self.column_size,
                                 cbar_ticks=self.cbar_ticks,
                                 density=self.density,
                                 density_cmap=self.density_cmap,
                                 bw=self.bw,
                                 density_alpha=self.density_alpha,
                                 dot_order=self.dot_order,
                                 inverse_dot_order=self.inverse_dot_order,
                                 order_cbar_labels=self.order_cbar_labels,
                                 vmin=self.vmin,
                                 vmax=self.vmax,
                                 grid=self.grid,
                                 markerscale=self.markerscale,
                                 labelsize=self.labelsize,
                                 bins=self.bins,
                                 axis_border=self.axis_border,
                                 showxticks=self.showxticks,
                                 showyticks=self.showyticks,
                                 legend_loc=self.legend_loc,
                                 color_dtype=self.color_dtype[count])
                    self.make_bottom_plot(ax_bottom=ax1, c_in=c)
                    row_count_1 += 1
                    row_count_2 += 1
                    if row_count_1 % self.column_size == 0:
                        row_count_1 += self.column_size
                        row_count_2 += self.column_size
            count += 1
        ## Set the extra plots to be invisible
        if self.clusters is not None:
            while row_count_1 % self.column_size != 0:
                ax0 = axs[row_count_1]
                ax1 = axs[row_count_2]
                ax0.set_visible(False)
                ax1.set_visible(False)
                row_count_1 += 1
                row_count_2 += 1
        if save_plot != False and self.show_cbar == True:
            fig.savefig(self.filename,
                        bbox_inches='tight',
                        dpi=self.dpi,
                        transparent=self.fig_transparent)
        elif save_plot != False and self.show_cbar == False:
            fig.savefig(self.filename,
                        bbox_inches='tight',
                        dpi=self.dpi,
                        transparent=self.fig_transparent)
        else:
            return (fig, axs)

    def make_bottom_plot(self, ax_bottom, c_in=None):
        '''
        We make a bottom violin plot based on cluster assignment.
        Calculate the median and quantiles of provided values based on
        provided groupings and make a violin plot.
        '''

        if c_in is None:
            c_in = self.c
        if len(np.unique(self.clusters)) == 1:
            color_by_clusters = [self.c]
            c = None
            labels = None
        else:
            color_by_clusters = []
            labels = []
            c = []
            if self.sort_by_str is not None:
                for u in self.sort_by_str:
                    curr_cluster_val = c_in[self.clusters == u]
                    color_by_clusters.append(curr_cluster_val)
                    c.append([np.mean(curr_cluster_val)])
                    labels.append(u)
            else:
                for u in np.unique(self.clusters):
                    curr_cluster_val = c_in[self.clusters == u]
                    color_by_clusters.append(curr_cluster_val)
                    c.append([np.mean(curr_cluster_val)])
                    labels.append(u)
            c = np.array(c)
            c = c / np.max(c)
        make_violin(X=color_by_clusters,
                    ax=ax_bottom,
                    labels=labels,
                    c=c,
                    s=self.s * 3,
                    fontsize=(self.fontsize - 2),
                    cmap=self.cmap,
                    rotation=self.rotation,
                    grid=False,
                    mean=True,
                    sina=self.sina,
                    return_obj=False,
                    vert=True)

class umap(scatter):
    def __init__(self,
                 X,
                 c=None,
                 filename='',
                 figsize=(3, 3),
                 cmap=None,
                 s=5,
                 title=True,
                 legend_title=None,
                 facecolor='white',
                 xlabel=None,
                 ylabel=None,
                 show_cbar=True,
                 alpha=1,
                 fontsize=13,
                 hspace=0.15,
                 readable_titles=False,
                 wspace=0.05,
                 insertion_spacer=0,
                 clusters=None,
                 rotation=0,
                 ncols=5,
                 sort_by_str=None,
                 sina=False,
                 run=True,
                 cbar_ticks=True,
                 density=False,
                 edgecolor=None,
                 linewidths=None,
                 density_cmap='Reds',
                 bw='scott',
                 density_alpha=0.6,
                 dot_order='ordered',
                 inverse_dot_order=False,
                 order_cbar_labels=[],
                 adj=None,
                 showxticks=False,
                 showyticks=False,
                 ylim=None,
                 xlim=None,
                 legend_loc=None,
                 vmin=None,
                 vmax=None,
                 stratify=False,
                 grid=False,
                 labelsize=12,
                 markerscale=None,
                 bins=300,
                 axis_border='off',
                 force_discrete_cbar=False,
                 sharex=True,
                 sharey=True,
                 fn='',
                 layer=None,
                 symbol=None,
                 components=[0, 1],
                 dpi=300,
                 **args):
        self.label_dtype = []
        if isinstance(X, anndata.AnnData):
            self._is_anndata = True
            if layer is not None:
                self.X = X.obsm['X_'+layer][:, components]
            else:
                if 'X_umap' in X.obsm.keys():
                    self.X = X.obsm['X_umap'][:, components]
                    layer = 'UMAP'
                elif 'X_pca' in X.obsm.keys():
                    self.X = X.obsm['X_pca'][:, components]
                    layer = 'PCA'
                else:
                    raise ValueError('No valid layer is provided')
            if isinstance(c, (str, list, pd.core.series.Series, np.ndarray, pd.core.indexes.base.Index)):
                if title:
                    try:
                        if symbol is not None:
                            title = X.var.loc[c, symbol]
                        elif stratify:
                            title = sorted(set(X.obs[c].values))
                        elif 'symbol' in X.var.columns:
                            title =  X.var.loc[c, 'symbol']
                        elif 'symbols' in X.var.columns:
                            title =  X.var.loc[c, 'symbols']
                        elif 'gene_symbol' in X.var.columns:
                            title =  X.var.loc[c, 'gene_symbol']
                        elif 'gene_symbols' in X.var.columns:
                            title =  X.var.loc[c, 'gene_symbols']
                        elif c is not None:
                            title = c
                        else:
                            title = ''
                    except:
                        title = ''
                if type(c) is str:
                    c = [c]
                _c = []
                for _plot_num, ci in enumerate(c):
                    assert ci in X.obs.columns or ci in X.var_names, '{:} is not a valid observation or feature key'.format(c)
                    if ci in X.obs.columns:
                        _cvals = X.obs[ci].values
                        _c.append(_cvals)
                        _c_dtype = X.obs[ci].values.dtype
                        self.label_dtype.append(_c_dtype)
                    else:
                        if sp.sparse.issparse(X.X):
                            _cvals = X[:, ci].X.A.flatten()
                        else:
                            _cvals = X[:, ci].X.flatten()
                        _c.append(_cvals)
                        _c_dtype = 'float'
                        self.label_dtype.append(_c_dtype)
                if len(c) <= 1:
                    c = np.array(_c[0])
                else:
                    c = np.array(_c)
            if type(adj) is str:
                adj = X.uns[adj]
            elif adj is True:
                if 'connectivities' in X.obsp.keys():
                    adj = X.obsp['connectivities']
                else:
                    adj = X.uns['neighbors']['connectivities']
            if type(clusters) is str:
                clusters = X.obs[clusters].values
        else:
            self._is_anndata = False
            self.X = X
        if c is None:
            self.c = np.zeros(np.shape(self.X)[0])
            show_cbar = False
        else:
            ndim_c = np.ndim(c)
            if adj is not None:
                assert np.shape(adj)[0] == np.shape(adj)[1],\
                        'Adjacency matrix should have syemmetric length'
                assert np.shape(adj)[0] == np.shape(self.X)[0],\
                        'Adjacency matrix should have the same length as\
                        number of samples'
                if ndim_c == 1:
                    self.c = np.array(adj.dot(c) / np.sum(adj, axis=1).reshape(
                        np.shape(c)[0], 1).T)[0]
                else:
                    self.c = average_by_nnm(c.T, adj).T
            else:
                if isinstance(c, pd.Categorical):
                    c = c.astype(str)
                self.c = np.array(c)
                if self._is_anndata is False:
                    for ci in self.c:
                        if type(ci) is str:
                            self.label_dtype.append('str')
                        else:
                            self.label_dtype.append(ci.dtype)
        self.clusters = clusters
        if self.clusters is not None:
            assert np.shape(clusters)[0] == np.shape(self.X)[0],\
                    'the length of clusters must equal the length of input\
                    samples'
            figsize = (figsize[0], figsize[1] + figsize[0] / 4)
            self.clusters = np.array(clusters)
        self.show_cbar = show_cbar
        if fn != '':
            self.filename = fn
        else:
            self.filename = filename
        self.legend_loc = legend_loc
        self.legend_title = legend_title
        if self.filename != '':
            self.save_plot = True
            dirname = os.path.dirname(os.path.abspath(self.filename))
            if os.path.isdir(dirname) == False:
                os.mkdir(dirname)
        else:
            self.save_plot = False
        self.column_size = int(ncols)
        self.figsize = figsize
        if stratify == True and cmap is None:
            self.cmap = {x:'C'+str(y) for y, x in enumerate(sorted(set(self.c)))}
        else:
            self.cmap = cmap
        if type(c) is not str and type(title) is bool:
            self.title = ''
        else:
            self.title = title
        self.facecolor = facecolor
        if xlabel or xlabel == '':
            self.xlabel = xlabel
        else:
            if layer:
                self.xlabel = layer.split('_')[-1].upper()+str(components[0]+1)
            else:
                self.xlabel = None
        if ylabel or xlabel == '':
            self.ylabel = ylabel
        else:
            if layer:
                self.ylabel = layer.upper().split('_')[-1]+str(components[1]+1)
            else:
                self.ylabel = None
        self.alpha = alpha
        self.fontsize = fontsize
        self.insertion_spacer = insertion_spacer
        self.grid = grid
        self.markerscale = markerscale
        self.labelsize = labelsize
        if self.insertion_spacer > 0:
            self.title = make_titles_readable(self.title,
                                    insertion_spacer=self.insertion_spacer)
        self.s = s
        self.hspace = hspace
        self.rotation = rotation
        self.sina = sina
        self.cbar_ticks = cbar_ticks
        self.density = density
        self.bw = bw
        self.density_alpha = density_alpha
        self.fig_transparent = False
        self.dot_order = dot_order
        self.inverse_dot_order = inverse_dot_order
        self.order_cbar_labels = order_cbar_labels
        self.edgecolor = edgecolor
        self.linewidths = linewidths
        self.showxticks = showxticks
        self.showyticks = showyticks
        self.bins = bins
        self.args = args
        self.ylim = ylim
        self.xlim = xlim
        self.vmin = vmin
        self.vmax = vmax
        self.stratify = stratify
        self.axis_border = axis_border
        self.force_discrete_cbar = force_discrete_cbar
        self.sharex = sharex
        self.sharey = sharey
        self.dpi = dpi
        if self.facecolor == 'w' or self.facecolor == 'white':
            self.fig_transparent = True
        if self.density == False and density_cmap == 'None':
            self.density_cmap = self.cmap
        else:
            self.density_cmap = density_cmap
        if sort_by_str is not None and self.clusters is not None:
            assert np.in1d(sort_by_str, np.unique(self.clusters)).sum() ==\
                    len(np.unique(self.clusters)), 'The provided sort_by_str\
                    must contain all unique cluster values'

            self.sort_by_str = np.array(sort_by_str)
        else:
            self.sort_by_str = sort_by_str
        if self.show_cbar == True:
            self.wspace = (wspace + 0.15)
        else:
            self.wspace = wspace
        self.readable_titles = readable_titles
        plt.rcParams["axes.edgecolor"] = "0.15"
        plt.rcParams["axes.linewidth"] = 1.25
        if run == True:
            self.run()


def scatter_plot(X=None,
                 c=None,
                 count=None,
                 xlabel='',
                 ylabel='',
                 fontsize=12,
                 labelsize=12,
                 title='',
                 number_of_plots=1,
                 facecolor='w',
                 cmap=None,
                 s=5,
                 alpha=1,
                 show_cbar=True,
                 column_size=1,
                 cbar_ticks=True,
                 density=False,
                 density_cmap='hot',
                 bw='scott',
                 density_alpha=0.6,
                 edgecolor=None,
                 linewidths=None,
                 dot_order='ordered',
                 inverse_dot_order=False,
                 order_cbar_labels=[],
                 showxticks=True,
                 showyticks=True,
                 textlabel=None,
                 connectivity=None,
                 hide_color=False,
                 legend_loc=None,
                 legend_title=None,
                 rasterized=None,
                 vmin=None,
                 vmax=None,
                 grid=False,
                 markerscale=None,
                 bins=300,
                 mode='scatter',
                 axis_border='on',
                 warning=True,
                 color_dtype=None,
                 zorder=None,
                 ax=None,
                 fig=None,
                 x=None,
                 y=None):
    if X is None:
        x = np.array(x)
        y = np.array(y)
        assert isinstance(x, list) or isinstance(x, np.ndarray), 'supplied x must be either of list or array type'
        assert isinstance(y, list) or isinstance(y, np.ndarray), 'supplied y must be either of list or array type'
        X = np.array([x, y]).T
    if count is None:
        count = 0
    if c is None:
        c = np.ones(np.shape(X)[0])
    if not color_dtype:
        color_dtype = c.dtype
    if showyticks == False:
        ax.set_yticks([])
    if showxticks == False:
        ax.set_xticks([])
    if hide_color == True:
        show_cbar = False
    if rasterized is None:
        if np.shape(X)[0] > 100:
            rasterized = True
        else:
            rasterized = False
    if axis_border!='sc':
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
    else:
        ax.set_xlabel(xlabel, fontsize=fontsize, x=0.1, ha='left')
        ax.set_ylabel(ylabel, fontsize=fontsize, y=0.1, ha='left')
    ax.set_title(title, fontsize=fontsize)
    ax.set_facecolor(facecolor)
    ax.tick_params(labelsize=labelsize)
    if grid is False:
        ax.grid(False)
    else:
        ax.grid(color='k', alpha=0.6, ls='--')
    if color_dtype in ['float', 'float32', 'float64', 'int32', 'int64']:
        if dot_order == 'ordered':
            if inverse_dot_order == True:
                sorted_id = np.argsort(-c)
            else:
                sorted_id = np.argsort(c)
        else:
            sorted_id = np.arange(len(c))
        if cmap is None:
            cmap = 'viridis'
        artist = ax.scatter(X[sorted_id, 0],
                            X[sorted_id, 1],
                            c=c[sorted_id],
                            cmap=cmap,
                            s=s,
                            alpha=alpha,
                            edgecolor=edgecolor,
                            linewidths=linewidths,
                            rasterized=rasterized,
                            zorder=zorder,
                            vmin=vmin,
                            vmax=vmax)
        if show_cbar == True:
            fmt = FormatScalarFormatter("%.1f")
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar = fig.colorbar(artist, cax=cax, format=fmt)
            if np.max(cbar.get_ticks()) > 100:
                cbar.formatter.set_powerlimits((0, 0))
            cbar.ax.yaxis.set_offset_position('left')
            cbar.ax.tick_params(labelsize=(fontsize - 2))
            if legend_title is not None:
                cbar.ax.set_title(legend_title, fontsize=(fontsize-2))
            cbar.update_ticks()
            if cbar_ticks == False:
                cbar.set_ticks([])
    else:
        labels = c
        unique_labels = np.unique(labels)
        if isinstance(dot_order, (list, np.ndarray)):
            unique_labels = np.array(dot_order)
        if len(unique_labels) > 200 and warning:
            raise ValueError('You have more than 200 unique classes in the color labels, are you sure these are discrete colors?')
        if inverse_dot_order == True:
            unique_labels = unique_labels[::-1]
        if cmap is None and len(unique_labels) <= 10:
            cmap = dict([(y, 'C'+str(x)) for x, y in enumerate(sorted(set(
                unique_labels)))])
        elif type(cmap) is not dict:
            tick_dictionary = dict([
                (y, x) for x, y in enumerate(sorted(set(unique_labels)))
            ])
            c = np.array([tick_dictionary[x] for x in np.unique(labels)])
            minima = min(c)
            maxima = max(c)
            norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
            mapper = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        if hide_color == False:
            if len(order_cbar_labels) > 0:
                for key in order_cbar_labels:
                    assert str(key) in labels, 'Provided sort label %s\
                    is not in the labels' % (str(key))
                    if type(cmap) is dict:
                        assert key in cmap.keys(), 'Provided sort label %s\
                        is not in the colormap dictionary' % (str(key))
                        ax.scatter(x=X[labels == key, 0],
                                   y=X[labels == key, 1],
                                   label=key,
                                   color=cmap[key],
                                   s=s,
                                   alpha=alpha,
                                   edgecolor=edgecolor,
                                   linewidths=linewidths,
                                   zorder=zorder,
                                   rasterized=rasterized)
                    else:
                        ax.scatter(x=X[labels == key, 0],
                                   y=X[labels == key, 1],
                                   label=key,
                                   color=mapper.to_rgba(tick_dictionary[key]),
                                   s=s,
                                   alpha=alpha,
                                   edgecolor=edgecolor,
                                   linewidths=linewidths,
                                   zorder=zorder,
                                   rasterized=rasterized)
            else:
                for key in unique_labels:
                    if type(cmap) is dict:
                        assert key in cmap.keys(), 'Provided sort label %s\
                        is not in the colormap dictionary' % (str(key))
                        ax.scatter(x=X[labels == key, 0],
                                   y=X[labels == key, 1],
                                   label=key,
                                   color=cmap[key],
                                   s=s,
                                   alpha=alpha,
                                   edgecolor=edgecolor,
                                   linewidths=linewidths,
                                   zorder=zorder,
                                   rasterized=rasterized)
                    else:
                        ax.scatter(x=X[labels == key, 0],
                                   y=X[labels == key, 1],
                                   label=key,
                                   color=mapper.to_rgba(tick_dictionary[key]),
                                   s=s,
                                   alpha=alpha,
                                   edgecolor=edgecolor,
                                   linewidths=linewidths,
                                   zorder=zorder,
                                   rasterized=rasterized)
        if not markerscale:
            markerscale = int(np.sqrt(10 / s) * 3)
        if connectivity is not None:
            lw = s / 10
            central_pos = make_mean_pos(X, clusters=labels)
            for ind, u in enumerate(central_pos.keys()):
                ax.scatter(central_pos[u][0], central_pos[u][1],\
                          s=s*5, cmap='Spectral', c='k')
            draw_edge(central_pos, connectivity=connectivity, ax=ax, max_lw=5)
        if show_cbar != False:
            if legend_loc is not None:
                if type(legend_loc) is str:
                    allowed_loc = [
                        'upper left', 'upper right', 'bottom left', 'bottom right'
                    ]
                    assert legend_loc in allowed_loc, 'legend_loc must be one of\
                    the available choices: {:}'.format(allowed_loc)
                    if legend_loc == 'upper left':
                        location = (0.05, 0.95)
                    elif legend_loc == 'upper right':
                        location = (0.95, 0.95)
                    elif legend_loc == 'bottom left':
                        location = (0.05, 0.05)
                    elif legend_loc == 'bottom right':
                        location = (0.95, 0.05)
                else:
                    location = legend_loc
                lgnd = ax.legend(fontsize=(fontsize - 2),
                                 loc=location,
                                 fancybox=False,
                                 shadow=False,
                                 edgecolor='black',
                                 frameon=False,
                                 markerscale=markerscale,
                                 borderaxespad=0.)
            else:
                lgnd = ax.legend(fontsize=(fontsize - 2),
                                 loc=2,
                                 bbox_to_anchor=(1.05, 1),
                                 fancybox=False,
                                 shadow=False,
                                 edgecolor='black',
                                 frameon=False,
                                 markerscale=markerscale,
                                 borderaxespad=0.)
            if legend_title is not None:
                lgnd.set_title(legend_title, prop={"size":(fontsize-2)})
    if axis_border != 'on':
        if axis_border == 'off':
            ax.axis('off')
        elif axis_border == 'sc':
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            xdist = (xmax - xmin)/5
            ydist = (ymax - ymin)/5
            ax.spines['bottom'].set_bounds(xmin, xmin+xdist)
            ax.spines['left'].set_bounds(ymin, ymin+ydist)

    if density == True:
        density2d(X, ax, bins=bins, xlim=ax.get_xlim(), ylim=ax.get_ylim(),
                s=s, cmap=density_cmap, mode=mode, mesh_order='top', alpha=alpha)


def make_titles_readable(titles, insertion_spacer=3):
    '''
    This is to shorten gene names and replace characters such as % and _
    It also inserts a linebreaker symbol at every insertion_spacer increment
    '''
    title_v = []
    for i in titles:
        title = str(i).replace('_', ' ').split(' ')
        for eid, i in enumerate(
                range(insertion_spacer - 1,
                      len(title) - 1, insertion_spacer)):
            title[i] += '\n'
        title = ' '.join(title)
        title_v.append(title.rstrip())
    return (np.array(title_v))


def set_axis_style(ax, labels, fontsize=12, rotation=0):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(0, len(labels)))
    ax.set_xticklabels(labels, fontsize=fontsize, rotation=rotation)


def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return (lower_adjacent_value, upper_adjacent_value)


def make_violin(X,
                figsize=(4, 4),
                ax=None,
                labels=None,
                c=None,
                fontsize=12,
                rotation=0,
                s=50,
                cmap='gnuplot',
                xlabel='',
                ylabel='',
                grid=False,
                mean=False,
                showyticks=False,
                title='',
                sina=False,
                sample_threshold=50,
                random_seed=True,
                facecolor='white',
                filename='',
                transparent=True,
                return_obj=True,
                vert=True,
                logy=False):
    X = np.array(X)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    number_of_samples = np.shape(X)[0]
    if c is None:
        if number_of_samples == 1:
            c = [1]
        else:
            c = np.ones(np.shape(X)[0])
    if number_of_samples == 1:
        inds = [0]
    else:
        inds = np.arange(0, np.shape(X)[0])
    if labels is None:
        labels = inds
    if title != '':
        ax.set_title(title, fontsize=fontsize, y=1.02)
    ax.tick_params(labelsize=fontsize)
    quartile1 = np.zeros(number_of_samples)
    medians = quartile1.copy()
    quartile3 = quartile1.copy()
    quant_median = [np.percentile(x, [25, 50, 75]) for x in X]
    mean_v = np.array([np.nanmean(x) for x in X])
    for eid, i in enumerate(quant_median):
        quartile1[eid] = i[0]
        medians[eid] = i[1]
        quartile3[eid] = i[2]
    whiskers = np.array([adjacent_values(sorted_array, q1, q3) for \
                         sorted_array, q1, q3 in zip(X, quartile1, quartile3)])
    whiskersMin, whiskersMax = (whiskers[:, 0], whiskers[:, 1])
    ax.set_facecolor(facecolor)
    ax.grid(False)
    if sina == False:
        if mean == True:
            ax.scatter(inds,
                       mean_v,
                       marker='o',
                       color='white',
                       edgecolor='k',
                       s=s,
                       zorder=3,
                       lw=1)
        else:
            ax.scatter(inds,
                       medians,
                       marker='o',
                       color='white',
                       edgecolor='k',
                       s=s,
                       zorder=3,
                       lw=1)
    else:
        x_inds = []
        swarmplot_y = []
        for i, j in zip(X, inds):
            nm = np.shape(i)[0]
            if nm > sample_threshold:
                swarmplot_y.append(
                    np.random.choice(i, sample_threshold, replace=False))
                x_inds.append(np.repeat(j, sample_threshold))
            else:
                swarmplot_y.append(i)
                x_inds.append(np.repeat(j, nm))
        swarmplot_y = np.concatenate(swarmplot_y)
        sns.swarmplot(x=np.concatenate(x_inds),
                      y=swarmplot_y,
                      ax=ax,
                      color='k',
                      s=int(s / 10),
                      alpha=0.7)
    parts = ax.violinplot([x for x in X],
                          positions=inds,
                          showmeans=False,
                          showmedians=False,
                          showextrema=False,
                          vert=vert)
    for pc, color in zip(parts['bodies'], c):
        pc.set_facecolor(plt.get_cmap(cmap)(color))
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
        pc.set_alpha(1)
    ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=3)
    ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)
    set_axis_style(ax, labels, fontsize=fontsize, rotation=rotation)
    ax.set_xlabel(xlabel, fontsize=fontsize, labelpad=5)
    ax.set_ylabel(ylabel, fontsize=fontsize, labelpad=5)
    if logy == True:
        ax.set_yscale('log')
        ymin = (np.max([np.max(x) for x in X]))
        ymax = (np.min([np.min(x) for x in X]))
        if ymax < 1:
            yticks = [x for x in ax.get_yticklabels()]
    if showyticks == False:
        ax.set_yticks([])
    if grid == True:
        ax.yaxis.grid(color='k', linestyle=':', lw=1)
    if filename != '':
        fig.savefig(filename, bbox_inches='tight', transparent=transparent)
    elif return_obj == True:
        return (fig, ax)


class stratify():
    def __init__(self,
                 X,
                 labels,
                 c=None,
                 binsize=10,
                 figsize=(3, 3),
                 filename='',
                 plot=True,
                 ax=None,
                 ylabel='',
                 xlabel='',
                 order_cbar_labels=[]):
        assert binsize % 2 == 0, 'binsize must be an even number'
        self.binsize = binsize
        self.X = X
        self.figsize = figsize
        self.labels = labels
        self.c = c
        self.filename = filename
        self.plot = plot
        self.ax = ax
        self.ylabel = ylabel
        self.xlabel = xlabel
        self.order_cbar_labels = order_cbar_labels
        self.bin_it()
        if self.plot == True:
            self.run()

    def bin_it(self, style='Percent'):
        length = np.shape(self.X)[0]
        df = pd.DataFrame(0,
                          index=np.unique(self.labels),
                          columns=range(length))
        for i in range(length):
            if i < int(
                    self.binsize / 2) or i > (length - int(self.binsize / 2)):
                if i < self.binsize:
                    l, n = 0, int(self.binsize / 2)
                else:
                    l, n = length - int(self.binsize / 2), length
            else:
                l, n = i - self.binsize, i + self.binsize
            lb, cnt = np.unique(self.labels[l:n], return_counts=True)
            df.loc[lb, i] = cnt
        self.df = df
        if style == 'Percent':
            self.df = self.df.div(self.df.sum(axis=0), axis=1) * 100
        elif style == 'Frequency':
            self.df = self.df.div(self.df.sum(axis=0), axis=1)
        if len(self.order_cbar_labels) > 0:
            self.df = self.df.loc[self.order_cbar_labels]

    def run(self):
        external_plot_control = True
        if self.ax is None:
            self.figure, self.ax = plt.subplots(figsize=self.figsize)
            external_plot_control = False
        else:
            self.ax = self.ax.twinx()
        x = np.arange(1, np.shape(self.X)[0] + 1)
        if external_plot_control == True:
            self.ax.yaxis.tick_right()
        for i in self.df.index.values:
            if self.c is None:
                self.ax.plot(x, self.df.loc[i].values, label=i)
            else:
                self.ax.plot(x,
                             self.df.loc[i].values,
                             label=i,
                             color=self.c[i])
        self.ax.set_ylabel(self.ylabel, rotation=0)
        if self.xlabel != '':
            self.ax.set_xlabel(self.xlabel)
        if self.filename != '' and external_plot_control == False:
            self.figure.savefig(self.filename,
                             bbox_inches='tight',
                             transparent=True)

def annotate(genenames, xs, ys, ax, gkey=None, adjust_text=True, fontsize=6,
        gstyle=1, color='k'):
    genenames = list(genenames)
    texts = []
    for i, g in enumerate(genenames):
        if gstyle==1:
            texts.append(ax.text(xs[i], ys[i], g, fontsize=fontsize, color=color))
        elif gstyle==2:
            texts.append(ax.text(xs[i], ys[i], g,
                fontsize=fontsize, bbox=dict(alpha=0.1, boxstyle='round', color=color)))
        else:
            raise ValueError('invalid gstyle choice')
    if adjust_text:
        fix_annotation(texts, ax=ax)

def text(x, y, s, ax, adjust_annotation=True, lw=1, **args):
    if isinstance(x, list) is False and isinstance(x, np.ndarray) is False:
        ax.text(x, y, s, **args)
    else:
        texts = []
        for i, xi, in enumerate(x):
            if isinstance(y, list) is False and isinstance(y, np.ndarray) is False:
                yi = y
            else:
                yi = y[i]
            if isinstance(s, list) is False and isinstance(s, np.ndarray) is False:
                si = s
            else:
                si = s[i]
            texts.append(ax.text(xi, yi, si, **args))
        if adjust_annotation:
            fix_annotation(texts, lw=lw)
        else:
            return texts

def fix_annotation(texts, lw=0.5, **args):
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=lw),
            **args)

def average_by_nnm(X, nnm):
    if sp.sparse.issparse(nnm):
        D_avg = (nnm.multiply(1/nnm.sum(1).A)).dot(X)
    else:
        D_avg = (nnm.multiply(1/nnm.sum(1))).dot(X)
    return D_avg

class volcano():
    def __init__(self, D="dataframe", lfc='log2-fold', pv='adjusted-p', lfc_thr=1, pv_thr=0.05, color=("xkcd:green", "xkcd:red"), valpha=1,
                geneid=None, annotate=None, gfont=6, figsize=(3, 3), r=300, ar=0, s=3, markerdot="o",
                sign_line=False, gstyle=1, fn='', figtype='png', axtickfontsize=9,
                axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, title=None, ax=None, rasterized=True):
        if ax is None:
            ax_supplied = False
            ax = plt
        else:
            ax_supplied = True
        self.d = D.copy()
        min_x = self.d[pv][self.d[pv] > 0].min()
        self.d.loc[self.d[pv] == 0, pv] = min_x
        self._x = r'$log_{2}$(Fold Change)'
        if pv == 'adjusted-p':
            self._y = r'$-log_{10}$(FDR)'
        else:
            self._y = r'$-log_{10}$(P-value)'
        self.color = color
        self.d.loc[(self.d[lfc] >= lfc_thr) & (self.d[pv] < pv_thr), 'color'] = color[0]  # upregulated
        self.d.loc[(self.d[lfc] <= -lfc_thr) & (self.d[pv] < pv_thr), 'color'] = color[1]  # downregulated
        self.d['color'].fillna('grey', inplace=True)  # intermediate
        self.d['logpv'] = -(np.log10(self.d[pv]))
        self.rasterized = rasterized
        if ax_supplied is False:
            plt.subplots(figsize=figsize)
        ax.scatter(self.d[lfc], self.d['logpv'], c=self.d['color'], alpha=valpha, s=s, marker=markerdot, rasterized=self.rasterized)
        if sign_line:
            ax.axhline(y=-np.log10(pv_thr), linestyle='--', color='#7d7d7d', linewidth=1)
            ax.axvline(x=lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
            ax.axvline(x=-lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
        if type(annotate) is int:
            upregulated = self.d.loc[(self.d[pv] < pv_thr) & (self.d[lfc] > 0)].sort_values(pv, ascending=True).index.values[:annotate]
            downregulated = self.d.loc[(self.d[pv] < pv_thr) & (self.d[lfc] < 0)].sort_values(pv, ascending=True).index.values[:annotate]
            annotate = np.concatenate([upregulated, downregulated])
        self.geneplot(self.d, geneid, lfc, lfc_thr, pv_thr, annotate, gfont, pv, gstyle, ax, ax_supplied)
        if axxlabel:
            self._x = axxlabel
        if axylabel:
            self._y = axylabel
        if title:
            if ax_supplied:
                ax.set_title(title, fontsize=axlabelfontsize)
            else:
                ax.title(title, fontsize=axlabelfontsize)
        general.axis_labels(self._x, self._y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(fn, r)

    def geneplot(self, d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle, ax, ax_supplied=False):
        if genenames is not None:
            genenames = tuple(genenames)
            if geneid is None:
                features = d.index.values
            else:
                features = d[geneid].values
            texts = []
            if genenames == "deg":
                for i in np.unique(features):
                    if (d.loc[features == i, lfc].iloc[0] >= lfc_thr and d.loc[features == i, pv].iloc[0] < pv_thr) or \
                            (d.loc[features == i, lfc].iloc[0] <= -lfc_thr and d.loc[features == i, pv].iloc[0] < pv_thr):
                        if gstyle==1:
                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0], i,
                                          fontsize=gfont))
                        elif gstyle==2:

                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0], i,
                                fontsize=gfont, bbox=dict(alpha=0.1, boxstyle='round')))
                        else:
                            print("Error: invalid gstyle choice")
                            sys.exit(1)
            elif isinstance(genenames, (tuple, list, np.ndarray)):
                for i in np.unique(features):
                    if i in genenames:
                        if gstyle==1:
                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0], i,
                                          fontsize=gfont))
                        elif gstyle==2:

                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0], i,
                                fontsize=gfont, bbox=dict(alpha=0.1, boxstyle='round')))
                        else:
                            print("Error: invalid gstyle choice")
                            sys.exit(1)
            elif type(genenames) is dict:
                for i in np.unique(features):
                    if i in genenames:
                        if gstyle==1:
                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0],
                                          genenames[i], fontsize=gfont))
                        elif gstyle == 2:
                            texts.append(ax.text(d.loc[features == i, lfc].iloc[0],
                                d.loc[features == i, 'logpv'].iloc[0], i,
                                fontsize=gfont, bbox=dict(alpha=0.1, boxstyle='round')))
                        else:
                            print("Error: invalid gstyle choice")
                            sys.exit(1)
            if ax_supplied is False:
                adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1, alpha=0.5, shrinkA=5))
            else:
                adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=1, alpha=0.5, shrinkA=5), ax=ax)

    def involcano(self, d="dataframe", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"),
                  valpha=1, geneid=None, genenames=None, gfont=8, figsize=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=True, figtype='png', axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None):
        _x = r'$log_{2}$(Fold Change)'
        _y = r'$-log_{10}$(P-value)'
        color = color
        d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = color[0]  # upregulated
        d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = color[1]  # downregulated
        d['color'].fillna('grey', inplace=True)  # intermediate
        d['logpv'] = -(np.log10(d[pv]))

        plt.subplots(figsize=figsize)
        plt.scatter(d[lfc], d['logpv'], c=d['color'], alpha=valpha, s=dotsize, marker=markerdot, rasterized=self.rasterized)
        gene_exp.geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)
        plt.gca().invert_yaxis()
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        if xlm:
            print('Error: xlm not compatible with involcano')
            sys.exit(1)
        if ylm:
            print('Error: ylm not compatible with involcano')
            sys.exit(1)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, 'involcano')

    def ma(self, df="dataframe", lfc="logFC", ct_count="value1", st_count="value2", lfc_thr=1, valpha=1, dotsize=8,
           markerdot="o", figsize=(6, 5), r=300, show=True, color=("green", "red"), ar=90, figtype='png', axtickfontsize=9,
           axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
           axylabel=None, xlm=None, ylm=None):
        _x, _y = 'A', 'M'
        df.loc[(df[lfc] >= lfc_thr), 'color'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr), 'color'] = color[1]  # downregulated
        df['color'].fillna('grey', inplace=True)  # intermediate
        df['A'] = np.log2((df[ct_count] + df[st_count]) / 2)
        plt.subplots(figsize=figsize)
        plt.scatter(df['A'], df[lfc], c=df['color'], alpha=valpha, s=dotsize, marker=markerdot)
        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, 'ma')

    def hmap(self, df="dataframe", cmap="seismic", scale=True, figsize=(4, 6), clus=True, zscore=None, xlabel=True,
             ylabel=True, tickfont=(10, 10), r=300, show=True, figtype='png'):
        fig, hm = plt.subplots(figsize=figsize)
        if clus:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=figsize)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, 'heatmap')
        else:
            hm = sns.heatmap(df, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
            plt.xticks(fontsize=tickfont[0])
            plt.yticks(fontsize=tickfont[1])
            general.get_figure(show, r, figtype, 'heatmap')

class general():
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def get_figure(fn, r=300):
        if fn != '':
            plt.savefig(fn, bbox_inches='tight', dpi=r)

    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)

    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

