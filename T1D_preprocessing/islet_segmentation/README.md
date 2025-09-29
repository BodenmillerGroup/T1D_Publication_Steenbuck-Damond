# Islet segmentation model

Code to train an islet segmentation model based on U-Net.


## Training data

Training data was obtained by performing a pilot experiment on 16 patients
(pancreas donors), with 5 antibody panels applied on 5 consecutive sections.

Islet segmentation on the pilot experiment images was performed using
[ilastik](https://www.ilastik.org) and [CellProfiler](https://cellprofiler.org/),
as described in the [IMC segmentation pipeline](https://bodenmillergroup.github.io/ImcSegmentationPipeline)
and in a previous publication ([Damond et al. Cell Metab 2019](https://doi.org/10.1016/j.cmet.2018.11.014)).

The obtained images and masks were processed using the `data_formatting.ipynb`
jupyter notebook. This notebook is only provided for documentation and the
processed images and masks can directly be downloaded from zenodo
***UPDATE: ADD LINKS TO ZENODO REPOSITORY***


## Model training

The islet model segmentation was trained on the processed training data using
the `islet_segmentation.ipynb`. Running this notebook requires the helper
functions in `helpers.py`.

The generated model can be directly downloaded from zenodo
***UPDATE: ADD LINKS TO ZENODO REPOSITORY***. This model is used for islet
segmentation in `../processing/02_IsletSegmentation.ipynb`.


## Installation

To run the notebooks in this subfolder, it is recommended to use the dedicated
[docker container](https://hub.docker.com/r/ndamond/t1d_processing) to run the
jupyter notebooks in the `processing` directory. This container can be pulled
from docker hub with `docker pull ndamond/t1d_processing:0.0.2`.

Alternately, a conda environment can be installed by following the steps below.  
This may not work on all operating systems, however. Packages versions may also
change, see https://bodenmillergroup.github.io/steinbock/latest/install-python/
(section "Package version conflict") for updated instructions.
```
conda create -n t1d_processing python=3.8
conda activate t1d_processing
pip install --upgrade tensorflow==2.8.0
pip install --upgrade -r requirements_deepcell.txt
pip install --no-deps deepcell==0.12.0
pip install --upgrade -r requirements.txt
pip install --no-deps steinbock[all]
pip install opencv-python-headless
conda install pytorch torchvision torchaudio cpuonly -c pytorch
pip install segmentation-models-pytorch
pip install -U albumentations --no-binary qudida,albumentations
conda install jupyterlab
```
The `requirements.txt` and `requirements_deepcell.txt` files can be found in
the `../ext` folder of this repository.
