# Processing

This folder contains jupyter notebooks to extract images from IMC acquisitions,
segment cell and islets, and export data for downstream analysis in R.  

Detailed documentation is included in the notebooks.


## Installation

To run the notebooks in this subfolder, it is recommended to use the dedicated
[docker container](https://hub.docker.com/repository/docker/nathanste/t1d_preprocess/general) to run the
jupyter notebooks in the `processing` directory. This container can be pulled
from docker hub with `docker pull nathanste/t1d_preprocess:latest`.

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

## How to run the Jupyter notebooks from within Docker?

Not familiar with Docker? Read here how to run the jupyter notebooks
from within the Docker container.

Pre-requisites:
- Install Docker
- Clone this repository
- Download the raw data from Zenodo (**ADD LINK**)

Alternatively, use VSCode or other IDEs with the required extensions.

### Run the Docker container

With this command you can run a Docker container interactively.

`docker run -it -p <local_PORT>:<container_PORT> -v </LOCAL_path/to/t1d_preprocessing/repo/>:</path_in/DOCKER/> -v </LOCAL_path/to_raw_files/>:</path_in/DOCKER> <DOCKER_IMAGE>`

Let's break the options down:
`-it`: Runs the Docker container interactively.  
`-v` : Lets you interact with LOCAL files in the DOCKER file system via so-called volumes.  
       Thereby, you can interact with the raw IMC data in the Docker container.  
       Make sure to provide the LOCAL paths to **(a)** the raw IMC data and **(b)** to this Github repository.  
`-p` : Specifying a port name. 

For example:  
`docker run -it -p 8080:8080 -v /home/T1D_preprocessing:/home/T1D_preprocessing -v /mnt/central_nas/processing/:/home/processing/ nathanste/t1d_preprocess:latest`

### Run the Jupyter Notebook

After running the above command you should be "within" the Docker container.  
Check by examining the file system in the container (e.g. run `ls`).  
From within the Docker container, you can run this command to start a jupyter notebook:  

`jupyter notebook </path_docker/to_notebook.ipynb> --ip 0.0.0.0 --port <SAME_PORT> --no-browser --allow-root`

**Important:** Check that the port number behind `--port` matches the port number in the `docker run` command.  

For example:  
`jupyter notebook /home/T1D_preprocessing/processing/01_Preprocessing.ipynb --ip 0.0.0.0 --port 8080 --no-browser --allow-root`

### Access the Notebook
Open [http://localhost:8080/](http://localhost:8080/) in your browser.
