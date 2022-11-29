# Path: src/variables/granuloma_lib.py
# Define shared variables and common libraries used in all scripts
# Author: Yuan Xue (xuesoso@gmail.com)
# Github: https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
# Updated on 28 Nov 2022
# License: MIT

#### Import libraries ####
import os, sys, sklearn, warnings, getpass, copy, json, string, anndata, gzip, itertools
from pathlib import Path
__homepath__ = str(Path.home())
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path+'/..')
import numpy as np
import pandas as pd
import scipy as sp
from scipy.cluster import hierarchy
import scipy.io
from samalg import SAM
from nheatmap import nhm
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_hex
from sklearn import linear_model
from sklearn.cluster import KMeans
from sklearn.metrics import r2_score
import seaborn as sns
sns.set()
from umap import UMAP
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scanpy as sc
import scanpy.external as sce
import scvelo as scv
import statsmodels.api as sm
from glob import glob
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#### Input/Output directory paths ####
infd = script_path+'/../../../data/h5ad/'
out = script_path+'/../../../output/'
fdm = out+'Figures/Manuscript/'

#### Create output directory for figures ####
os.makedirs(fdm, exist_ok=True)

#### Figure color settings ####
cmap = {'A':'#787CCC', 'B':'#CC78A6', 'C':'#78cc9e', 'T':'#ccc878'}
dt = {'A':'WT', 'B':r'$\Delta$steE', 'C':'Control Ab', 'T':r'Anti-TNF'}
dtcmap = {dt[x]:cmap[x] for x in cmap}

#### Figure settings ####
output_formats = [".png", ".svg", ".pdf"]
savefig_args = {"dpi": 500, "bbox_inches": "tight", "pad_inches": 0.05, "transparent":True}

#### Editable textbox ####
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

#### A function to save figures in pdf, svg, and png ####
def sfig(fig, name, rasterized=False):
    for output_format in output_formats:
        fig.savefig(name + output_format, **savefig_args, rasterized=rasterized)
    return None
