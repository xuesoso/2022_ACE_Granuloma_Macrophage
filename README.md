# 2022_ACE_Granuloma_Macrophage

## Background

**This is a repository for the scRNA-seq analysis generated as part of the manuscript:**

Single-cell profiling identifies ACE+ granuloma macrophages as a non-permissive niche for intracellular bacteria during persistent Salmonella infection (2022). *Trung H. M. Pham†‡, Yuan Xue†, Susan M. Brewer, Kenneth E. Bernstein, Stephen R. Quake‡, Denise Monack‡.*

Legends:
†: co-first authors.
‡: co-corresponding authors.


## How to retrieve the datasets

**Processed datasets:**

- [Google drive deposit](https://drive.google.com/drive/folders/1ohx-A5gmWS42yG77ee6h7KLyaY4CINQV?usp=sharing)

- [GEO repository (GSE215880)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215880)

**Raw fastq:**

- [GEO repository (GSE215880)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215880)

**Notes:**
- The easiest way to access the processed data is through the [Google drive](https://drive.google.com/drive/folders/1ohx-A5gmWS42yG77ee6h7KLyaY4CINQV?usp=sharing)
.
- Download `data.zip` and unzip the contents to the data directory inside this repository:

    ```bash
    cd 2022_ACE_Granuloma_Macrophage/
    ```
    ```bash
    unzip data.zip -d ./data/
    ```

## Reproducing analysis in this repository

**1. Clone this repository to your local computer**

```bash
git clone https://github.com/xuesoso/2022_ACE_Granuloma_Macrophage
```



**2. Build the packages required by the code in this notebook**

- You have two options here:

    I. The easiest way is to build the Docker image with the Dockerfile provided in this repository and run the Jupyter notebook inside a container.

    - **Required: [Install Docker](https://docs.docker.com/get-docker/)**

    - Once you have Docker installed, navigate to the local directory of this Github repository:

    ```bash
    cd 2022_ACE_Granuloma_Macrophage/
    ```

    - Execute the script to build the image with pre-specified configurations:

    ```bash
    bash ./Docker/build_docker.sh
    ```

    - Execute the script to run this notebook under the Docker container:

    ```bash
    bash ./Docker/run_docker.sh
    ```

    - Navigate to the local address of the Jupyter notebook on your favorite browser, the default port passed is set to `8887`:

    ```bash
    firefox http://localhost:8887/notebooks/notebook/analysis_notebook.ipynb
    ```


    II. The second approach is to manually install the exact library versions. **This is not recommended as it involves navigating a complicated dependency graph**:

    - Install [Anaconda v4.8.3 with python 3.8](https://repo.anaconda.com/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh) on Linux system (tested on Ubuntu 18.04/Fedora 36).

    - Execute in shell: `conda install -c conda-forge install python=3.7.0`

    - Install the exact version of python packages with `pip` as documented in `requirements.txt`

    - *Note: You may have to install different dependent packages than the latest versions recommended by your platform's package manager.*

           
## Additioal notes
-------

**Supported operating system platforms:**
- MacOS (tested on Monterey v 12.5.1)
- Linux (tested on Fedora 36)

<details>
<summary>What is in each annotated data object</summary>

| data objects   |  descriptions |
| :---       |    :---   |
| sam_full.210505.h5ad | Processed 10X chromium v3.1 scRNA-seq data of splenocyte isolate. Contains samples collected from mouse infected by WT STm and dSTeE STm. Shown in Figure 1A-C. Shown in Figure S2A-B. |
| sam_full_velocyto.210505.h5ad | Contains same set of cells in "sam_full.210505.h5ad". Dataset is pre-processed with velocyto to yield RNA-velocity prediction. |
| sam_myeloidSubset.210505.h5ad | Myeloid sub-population derived from "sam_full.210505.h5ad". |
| sam_monocyteSubset.210505.h5ad | Monocyte macrophage sub-population derived from "sam_myeloidSubset.210505.h5ad". Shown in Figure 1D-G, Figure 2A-D, Figure S3A-C. |
| sam_monocyte_velocyto.210505.h5ad | Monocyte macrophage sub-population, same set of cells as in "sam_monocyteSubset.210505.h5ad". Dataset is pre-processed with velocyto to yield RNA-velocity prediction. Shown in Figure 2G. |
| sam_macrophageSubset_AB.210919.h5ad | Macrophage sub-population derived from "sam_monocyteSubset.210505.h5ad". Shown in Figure 3A-H. |
| harmony.sam_T_C_treatment.200119.sam_full.210505.h5ad | Processed 10X chromium v3.1 scRNA-seq data of splenocyte isolate. Contains samples collected from mouse treated with isotype control antibody or anti-TNFa antibody and then infected by WT STm. Data was aligned to "sam_full.210505.h5ad" using python implentation of harmony method (https://github.com/slowkow/harmonypy). |
| ABCT_SAM_momac.220828_review.h5ad | Derived from "sam_macrophageSubset_AB.210919.h5ad" and "harmony.sam_T_C_treatment.200119.sam_full.210505.h5ad". Contains samples collected from mouse infected by WT STm and dSTeE STm and from mouse treated with control antibody or anti-TNFa antibody. Both datasets have been subset to contain only the macrophage sub-population. "fig3_cell_type" and "fig3_leiden" have been updated to show cluster labels that are consistent with "sam_macrophageSubset_AB.210505.h5ad", which are the labels shown in Figure 6. "cell_type" and "leiden" reflect the original cluster labels generated by the analysis. |
| 201221_10X_velocyto_all.h5ad | This h5ad includes raw count values and velocyto estimate for all samples, including cells isolated from 4 x WT STm infected mice ("A"), 4 x dSTeE STm infected mice ("B"), 2 x isotype control antibody treated and WT STm infected mice ("C"), and 2 x anti-TNFa antibody treated and WT STm infected mice ("T"). |

</details>


