# Antibody titration

This folder contains a jupyter notebook that facilitates antibody titration for
Imaging Mass Cytometry experiments.  

To titrate antibodies, an antibody mix with starting concentrations was
pipetted and a serial dilution was performed. Here, the relative concentrations
used were 1500, 1000, 500, 250, and 125.

Tissue sections were then stained with these dilutions (one slide per dilution)
and imaged by imaging mass cytometry.

Raw IMC data were then processed using [steinbock](https://github.com/BodenmillerGroup/steinbock)
and exported to the [AnnData](https://anndata.readthedocs.io/en/latest/) format
for titration analysis using the `titration.ipynb` notebook provided in this
repository.

This notebook can be run using the files provided in the associated zenodo
repository ***ADD LINK***. Detailed instructions are provided within the
notebook.


## Installation

To run this script, a anaconda environment containing all required packages
can be installed using the `environment.yml` file.

The environment can be created using the following commands:
- `conda env create -f environment.yml` # Adapt the path to the environment.yml file if needed.  
- `conda activate titration`  
- `jupyter notebook`
