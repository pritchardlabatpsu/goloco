# goloco

goloco is a bioinformatics web application designed to enable genome-wide CRISPR loss-of-function predictions with small scale experiments of 100-1000 sgRNA subsets powered by lossy compression models. 

## What is lossy compression (loco)?
Trained with robust compendia of genome-wide CRISPR knockout experiments, our machine learning models identify cross predictive features that capture functional relationships between related genes and they self organize into growth regulatory networks. These networks can be compressed to only 100-1000 highly cross predictive nodes that predict the remaining ~17000-17900 genes with tunable loss of information, demonstrating the feasibility of compressing sgRNA libraries for genome-scale experiments. These lossy subsets make previously insurmountable experiments possible by generating genome-scale portraits of growth regulation with tiny pools of hundreds of sgRNAs. We hope to make functional genomic inquiry more effecient and scalable with this web application. Check out our previous [publication](https://www.nature.com/articles/s41467-022-28045-w) in Nature Communications describing our lossy compression algorithms in more detail [[1]](#1).

## Options for usage:
There are several options for using this application:
1. [Use the public version of goloco](http://goloco.herokuapp.com/) (**Recommended Overall**)
2. [Run goloco locally using Docker](#run-with-docker) (**Recommended Local**)
3. [Run goloco locally with Python](#run-with-python)

## Public Version:
Using the public version of this application is recommended for all users and is specifically designed to make this tool broadly accessible to those with limited or no prior coding knowledge. The application is hosted on [http://goloco.herokuapp.com](http://goloco.herokuapp.com/).

## Local Development:
Running this application on your local machine can be benefical to overcome CPU, memory, and storage limitations on the public server that limit genome-wide prediction runtimes and application responsiveness. If using either Docker or Python Dash to run this application locally, it is important to visit this manuals [data](#file-structure-and-data) section for instructions on how to download the core predictive models, store them locally, and modify the configuration variables to your local drives. This data step can be done immediately after cloning this repository for local development:

```bash
git clone https://github.com/pritchardlabatpsu/goloco.git
```

## Run with Docker:
If you wish to run goloco on your local machine, it is **recommended** to launch it as a Docker container which is fully configured to develop the application environment and launch the application services with little manual input. This procedure requires prior installation of [Docker desktop](https://www.docker.com/products/docker-desktop/).

### Step 1: Install Redis for Docker:
[Redis](https://redis.io/docs/getting-started/installation/) is an in-memory key-value store database used in this application as a message broker for long callback requests, native to Python Dash, with [celery](https://docs.celeryq.dev/en/stable/userguide/workers.html) backend workers processing requests. 

Open Docker desktop, search and pull the image for 'redis:latest', or run the following command in a terminal:

```bash
docker pull redis
```

### Step 2: Run Docker Compose:
Prior to launching this application, please see the [data](#file-structure-and-data) section for instructions on downloading the core predictive models and modifying configuration variables in the source code.

Clone this repository, navigate to the goloco directory that contains the Dockerfile and docker-compose.yml files, and run the following command in your terminal:

```bash
docker compose up
```

This previous step may take a while.

### Step 3: Open goloco:
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
conda create --name goloco -c conda-forge python=3.8
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
goloco uses a modification of the chronos software by the Broad Institute [[2]](#2) to generate psuedo-chronos gene effect scores from smaller scale CRISPR knock out experiments (i.e, lossy subset experiments). Visit the [chronos](https://github.com/broadinstitute/chronos) documentation for any specific installation instructions that may be relevant to your machine.

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
Redis is an in-memory key-value store database used in this application as a message broker for long callback requests native to Python Dash with celery backend workers processing requests. Visit the [Redis](https://redis.io/docs/getting-started/installation/) documentation for any specific installation instructions that may be relevant to your machine.

Run the following commands to install Redis for Linux distros or WSL2:

```bash
sudo apt-get update
sudo apt-get install redis-server
sudo systemctl enable redis-server.service
```

### Step 4. Launch goloco:
Prior to launching this application, please see the [data](#file-structure-and-data) section for instructions on downloading the core predictive models and modifying configuration variables in the source code.

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
celery -A app.celery_worker worker --loglevel=INFO
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

## File Structure and Data:
Some of the data for this application is already included in the data folder for this repository however the main gene models for making genome predictions and the predetermined clustering models need to be downloaded for local development. The file structure for this project is as follows:

```
.
├── assets             # empty folder for static images for web application.. no images in current deployment
├── cache              # empty folder for dynamic storage used by various callbacks and functions, i.e. temporary files generated during prediction
├── data               # folder for static data files used by application during various callbacks and functions
│   ├── df_crispr_19q4.pkl        # main dataframe from 
├── pages                  # contains source code for layouts and callbacks for each unique page hosted on the application
├── requirements.txt       # Python DEV USE: package specifications for pip installing requirements.. execute with "pip install -r requirements.txt"
├── app.py                 # Python DEV USE: main executable code with majority of dynamic functions, callbacks, and server initiation.. execute with "python app.py"
├── app_session.py         # code containing main class and asssociated functions for performing genome-wide prediction
├── celery_worker.py      # Python DEV USE: executable code to iniatite background celery workers.. execute with "celery -A celery_worker worker"
├── Dockerfile             # DOCKER DEV USE: executable code for creating single Docker container.. execute with "docker run"
├── docker-compose.yml     # DOCKER DEV USE: executable code using Dockerfile to create multicontainer Docker app.. execute with "docker compose up"
├── heroku.yml             # FOR PRODUCTION USE: heroku executable code for Heroku git deploy with Dockerfile to create multicontainer Heroku app
├── Procfile               # FOR PRODUCTION USE: heroku executable code for Heroku git deploy
├── .dockerignore          # instructs docker to ignore python cache files and env files
├── .gitignore             # instructs git to ignore python cache files and env files
```



## References

[1] Zhao B, Rao Y, Leighow S, O'Brien EP, Gilbert L, Pritchard JR. A pan-CRISPR analysis of mammalian cell specificity identifies ultra-compact sgRNA subsets for genome-scale experiments. Nat Commun. 2022 Feb 2;13(1):625. doi: 10.1038/s41467-022-28045-w. PMID: 35110534; PMCID: PMC8810922.

[2] Dempster, J.M., Boyle, I., Vazquez, F. et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol 22, 343 (2021). https://doi.org/10.1186/s13059-021-02540-7


