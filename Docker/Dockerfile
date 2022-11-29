FROM mambaorg/micromamba:0.8.2
ARG USER_ID
ARG GROUP_ID

#### Install packages with micromamba (much faster than anaconda) ####
RUN micromamba install -c plotly -c conda-forge \ 
                            python=3.7.0 \
                            pip=22.3.1 \
                            jupyter=1.0.0 \ 
                            plotly=5.11.0 \ 
                            numpy=1.21.6 \ 
                            scipy=1.7.3 \ 
                            pandas=1.1.5 \
                            umap-learn=0.3.10

RUN pip install adjustText==0.7.3 \
                anndata==0.8.0 \
                matplotlib==3.1.3 \
                DensityPlot==0.1.8 \
                scikit-learn==1.0.2 \
                statsmodels==0.13.5 \
                scanpy==1.6.0 \
                scvelo==0.1.25 \
                sam-algorithm==0.7.1 \
                seaborn==0.12.1 \
                Jinja2==3.1.2 \
                python-igraph==0.8.0

RUN micromamba install -c "conda-forge/label/cf202003" numba=0.48.0

#### Configure permission settings ####
RUN chmod ugo+rwx /root
RUN mkdir /jupyter && chmod ugo+rwx /jupyter
WORKDIR /jupyter/

#### Create a user with the same UID as the host user ####
RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
USER user

#### Run Jupyter Notebook backend ####
CMD jupyter notebook --port=8887 --no-browser --ip=0.0.0.0 --allow-root --NotebookApp.password="" --NotebookApp.token=""
