# goloco

goloco is a bioinformatics web application designed to enable genome-wide CRISPR loss-of-function predictions with small scale experiments of 100-1000 sgRNA subsets powered by lossy compression models. 

## What is lossy compression (loco)?
Trained with robust compendia of genome-wide CRISPR knockout experiments, our machine learning models identify cross predictive features that capture functional relationships between related genes and they self organize into growth regulatory networks. These networks can be compressed to only 100-1000 highly cross predictive nodes that predict the remaining ~17000-17900 genes with tunable loss of information, demonstrating the feasibility of compressing sgRNA libraries for genome-scale experiments. These lossy subsets make previously insurmountable experiments possible by generating genome-scale portraits of growth regulation with tiny pools of hundreds of sgRNAs. We hope to make functional genomic inquiry more efficient and scalable with this web application. Check out our previous [publication](https://www.nature.com/articles/s41467-022-28045-w) in Nature Communications describing our lossy compression algorithms in more detail [[1]](#1).

## Options for usage:
There are several options for using this application:
1. [Use the public version of goloco](http://goloco.herokuapp.com/) (**Recommended Overall**)
2. [Run goloco locally using Docker](#run-with-docker) (**Recommended Local**)
3. [Run goloco locally with Python](#run-with-python)

## Public Version:
Using the public version of this application is recommended for all users and is specifically designed to make this tool broadly accessible to those with limited or no prior coding knowledge. The application is hosted on [http://goloco.herokuapp.com](http://goloco.herokuapp.com/).

## Local Development:
Running this application on your local machine can be benefical to overcome CPU, memory, and storage limitations on the public server that limit genome-wide prediction runtimes and application responsiveness. If using either Docker or Python Dash to run this application locally, it is important to follow the instructions in the [data](#data) subsection to ensure models are accessible to the local app. First clone this repository:

```bash
git clone https://github.com/pritchardlabatpsu/goloco.git
```

### File Structure:
The cloned repository will have the following file structure:

```bash
.
├── assets                 # contains static images on web application
    ├──                            # images
├── cache                  # for temporary dynamic storage functions
    ├──                            # empty
├── ceres-infer            # empty directory which will serve as mount point for inference models
    ├──                            # empty
├── chronos                # chronos v2.0.5 from depmap git repo for direct installation in container
    ├──                            # installation files
├── data                   # contains static data files used by application, i.e. depmap19q4 release
    ├──                            # details described in data section
├── pages                  # contains source code for each unique webpage
    ├──                            # home.py, predict.py, overview.py, genes.py, zscore.py, clusters.py
├── requirements.txt       # PYTHON DEV USE: package specifications.. execute with "pip install -r requirements.txt" in venv
├── app.py                 # PYTHON DEV USE: main executable application server code.. execute with "python app.py or gunicorn"
├── celery_worker.py       # PYTHON DEV USE: executable background celery workers.. execute with "celery -A celery_worker worker"
├── app_session.py         # main class and asssociated functions for performing genome-wide prediction
├── Dockerfile             # DOCKER DEV USE: executable to create single Docker image.. execute with "docker run"
├── compose.dev.yml        # DOCKER DEV USE: executable to create multicontainer Docker app from scratch.. execute with "docker compose -f compose.dev.yml up"
├── compose.prod.yml       # DOCKER USE: executable to create multicontainer Docker app from public image.. execute with "docker compose -f compose.prod.yml up"
├── heroku.yml             # FOR PRODUCTION USE: heroku executable for git deploy to create multicontainer Docker based Heroku app
├── Procfile               # FOR PRODUCTION USE: heroku executable code for Heroku git deploy
├── .dockerignore          # instructs docker to ignore python cache files and env files
├── .gitignore             # instructs git to ignore python cache files and env files
├── .env                   # environmental variable required to run app in container
```

This repository is also publically available through Zenodo ([https://doi.org/10.5281/zenodo.14249073](https://doi.org/10.5281/zenodo.14249073)). This DOI does not contain the inference dataset as it is linked to the tagged releases from this Github Repo. See data section below for Zenodo library with dataset.

## Data:
Some data files are already included in the data folder of this repository. The application mainly uses [DepMap19Q4](https://depmap.org/portal/download/all/?releasename=DepMap+Public+19Q4) release files; our models were validated on this release [[1]](#1). Cell lines, CERES gene-effect scores, and definitions for essential, nonessential, and conditional genes are all consitent with DepMap19q4 [[2]](#2). 

Network files included in this repository are consistent with louvain communities determined in our previous work [[1]](#1).

### Download Models:
Genome-wide inference models and predetermined subclustering models of louvain communities need to be downloaded for local development and use:

The ceres-infer.zip file in the Zenodo library ([https://zenodo.org/records/14251390](https://zenodo.org/records/14251390)) contains the inference models. We recommend unzipping this file into a local directory outside of this repository as it will later be mounted onto a container. Unzipping this download somewhere in your local machine will create a new folder with the following subdirectories:

```bash
.
.
├── ceres-infer
    ├── L200_models      # contains all feature and pickled model files for each gene model
        ├──              
    ├── cluster_models   # contains all pickled clustering models for louvain communities
        ├── 
    ├── gw_predictions   # empty directory which will store genome wide predictions
        ├──              
    ├── gene_effects     # empty directory which will store gene effect score conversions
        ├── 
```

### Modify .env file:
The .env file will contain the following environmental variables specifying the PATH to the inference models from above and specifying the CPU type. Change the CERES_INFER_MODELS_PATH to the ceres-infer directory on your local machine prior to building the container. The directory will be mounted on the container during the build process. Change the CPU variable to either "intel" or "apple_silicon" depending on the CPU on your device.

```
CERES_INFER_MODELS_PATH=C:\Users\shasa\Desktop\ceres-infer <- change PATH to your local directory with inference models 
CPU=intel <- change to apple_silicon if CPU is apple based
```

## Run with Docker:
If you wish to run goloco on your local machine, it is **recommended** to launch it as a Docker container which is fully configured to develop the application environment and launch the application services with little manual input. This procedure requires prior installation of [Docker desktop](https://www.docker.com/products/docker-desktop/). Public prebuilt images are included in the Dockerhub under the name, "shossainova1\goloco-webapp". Images are available for intel CPUs tagged "shossainova1\goloco-webapp:intel" and apple silicon CPUs tagged "shossainova1\goloco-webapp:apple_silicon". Although not required prior to running step 1 below, this image can be pulled with the command ```docker pull shossainova1\goloco-webapp:{CPU}```, replacing {CPU} with the correct CPU tag.

### Step 1: Run Docker Compose:
Prior to launching this application, please see the [data](#data) section for instructions on downloading the core predictive models and modifying environmental variables in the source code.

Clone this repository, navigate to the goloco directory that contains the Dockerfile and compose.prod.yml files, and run the following command in your terminal:

```bash
docker compose -f compose.prod.yml up
```

This previous step may take a while. It will pull both the Redis image and the public webapp image from Dockerhub. It will also start all the app services and mount the inference models in your local directory to the container. [Redis](https://redis.io/docs/getting-started/installation/) is an in-memory key-value store database used in this application as a message broker for long callback requests, native to Python Dash, with [celery](https://docs.celeryq.dev/en/stable/userguide/workers.html) backend workers processing requests. 

**ALTERNATIVELY:** The container can be built from scratch on your local machine with the following command. This can be useful for developers who wish to change the app and rebuild the container. This can be done by running the following command:

```bash
docker compose -f compose.dev.yml up
```

### Step 2: Open goloco:
Once the previous step is completed, goloco should now be operational and serviced over port 8080 on your local machine. Open any web browser and navigate to the following address to access this application:

```bash
127.0.0.1:8080
```
or
```bash
localhost:8080
```

Your Docker desktop application will now have a fully contained version of the goloco app, which can be stopped and restarted at any time.


## Run with Python:
goloco can be installed and launched as a Python Dash application on your local machine. This procedure requires prior installation of [Anaconda](https://www.anaconda.com/) or [miniconda3](https://docs.conda.io/en/latest/miniconda.html).

### Step 1. Create Environment:
Clone this repository, navigate to the goloco main directory, and run the following commands to create and activate a new python virtual environment:

```bash
conda create --name goloco -c conda-forge python=3.9
conda activate goloco
```

This command can be used anytime to deactivate the virtual environment, however it is not intended for use prior to following the remaining steps:

```bash
conda deactivate
```

Run the following commands to install the basic required packages. Note that scikit-learn is installed manually to avoid conflicting dependencies specified by other packages:

```bash
pip install -r requirements.txt
pip install scikit-learn==0.22.1
```

### Step 2. Install CHRONOS:
goloco uses a modification of the chronos software by the Broad Institute [[3]](#3) to generate psuedo-chronos gene effect scores from smaller scale CRISPR knock out experiments (i.e, lossy subset experiments). Visit the [chronos](https://github.com/broadinstitute/chronos) documentation for any specific installation instructions that may be relevant to your machine.

In general, manual installation of chronos is required by running the following commands in the goloco main directory with the activated virtual environment:

```bash
git clone https://github.com/broadinstitute/chronos.git
```

Then navigate to the cloned chronos directory and run the following command:

```bash
python setup.py install
```
or
```bash
pip install .
```

### Step 3. Install Redis:
Redis is an in-memory key-value store database used in this application as a message broker for long callback requests, native to Python Dash, with celery backend workers processing requests. Visit the [Redis](https://redis.io/docs/getting-started/installation/) documentation for any specific installation instructions that may be relevant to your machine.

Run the following commands to install Redis for Linux distros or WSL2:

```bash
sudo apt-get update
sudo apt-get install redis-server
sudo systemctl enable redis-server.service
```

### Step 4. Launch goloco:
Prior to launching this application, please see the [data](#data) section for instructions on downloading the core predictive models and modifying configuration variables in the source code.

This app uses three main services to operate: Redis (message broker), Celery (backend worker), and Python Dash (main app). To launch these services open three separate terminals in the main goloco directory, activate the virtual environment in each, and run the following commands: 

#### Terminal 1. Start the Redis Service:
```bash
sudo service redis-server start
redis-cli
# you should see 127.0.0.1:6379>, then type 'monitor':
monitor
```

#### Terminal 2. Start the Celery Backend Service:
```bash
celery -A celery_worker worker --loglevel=INFO
```

#### Terminal 3. Start the Python Dash App Service:
```bash
python app.py
```

### Step 5. Open goloco:
Once the previous step is completed, goloco should now be operational and serviced over port 8080 on your local machine. Note that the port can be changed in the app.py source code if necessary or desired. Open any web browser and navigate to the following address to access this application:

```bash
127.0.0.1:8080
```
or
```bash
localhost:8080
```
The application can be stopped by interrupting or closing the processes in each of the service terminals and can be restarted by following Step 4.

## References

<a id="1">[1]</a> 
Zhao B, Rao Y, Leighow S, O'Brien EP, Gilbert L, Pritchard JR. A pan-CRISPR analysis of mammalian cell specificity identifies ultra-compact sgRNA subsets for genome-scale experiments. Nat Commun. 2022 Feb 2;13(1):625. doi: 10.1038/s41467-022-28045-w. PMID: 35110534; PMCID: PMC8810922.

<a id="2">[2]</a> 
Meyers RM, Bryan JG, McFarland JM, Weir BA, Sizemore AE, Xu H, Dharia NV, Montgomery PG, Cowley GS, Pantel S, Goodale A, Lee Y, Ali LD, Jiang G, Lubonja R, Harrington WF, Strickland M, Wu T, Hawes DC, Zhivich VA, Wyatt MR, Kalani Z, Chang JJ, Okamoto M, Stegmaier K, Golub TR, Boehm JS, Vazquez F, Root DE, Hahn WC, Tsherniak A. Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. Nat Genet. 2017 Dec;49(12):1779-1784. doi: 10.1038/ng.3984. Epub 2017 Oct 30. PMID: 29083409; PMCID: PMC5709193.

<a id="3">[3]</a> 
Dempster, J.M., Boyle, I., Vazquez, F. et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol 22, 343 (2021). https://doi.org/10.1186/s13059-021-02540-7

