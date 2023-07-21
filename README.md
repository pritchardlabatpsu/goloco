# goloco

goloco is a bioinformatics web application designed to enable genome-wide CRISPR loss-of-function predictions with small scale experiments of 100-1000 sgRNA subsets. Our machine learning models, trained with robust compendia of genome-wide CRISPR knockout experiments, identify predictive features that capture functional relationships between related genes and they self organize into growth regulatory networks. [[1]](#1). Information leveraged by highly cross predictive nodes, i.e. lossy compressed subsets of 100-1000 genes, in these networks can make previously insurmountable experiments possible by generating genome-scale portraits of growth regulation captured with tiny pools of hundreds of sgRNAs. With this tool we hope to make functional genomic inquiry more effiecient and scalable.

## Options for usage:
There are several options for using this application:
1. [Use the public version of goloco](http://goloco.herokuapp.com/) (**Recommended Overall**)
2. [Run goloco locally using Docker](#run-with-docker) (**Recommended Local**)
3. [Run goloco locally with Python](#run-with-python)

## Run with Docker:
Running this application on your local machine can be benefical to overcome CPU, memory, and storage limitations on the public server that pose limits to prediction runtimes and application responsiveness. If you wish to run goloco on your local machine, it is recommended to launch it as a Docker container which is fully configured to develop the application environment and launch the application services with little manual input. This procedure requires prior installation of [Docker desktop](https://www.docker.com/products/docker-desktop/).

### Step 1: Install Redis for Docker:
Open Docker desktop, search and pull the image for 'redis:latest', or run the following command in a terminal:

```bash
docker pull redis
```

### Step 2: Run Docker Compose:
Clone this repository, navigate to the goloco directory that contains the Dockerfile and docker-compose.yml files, and run the following command in your terminal:

```bash
docker compose up
```

This previous step may take a while.

### Step 3: Open goloco:
Once the previous step is completed, goloco should now be operational and serviced over port 8080 on your local machine. Open any web browser and type '127.0.0.1:8080' or 'localhost:8080' to access the application. Your Docker desktop application will now have a fully encolsed container for the goloco app which can be stopped and restarted and any time.


## Run with Python:
goloco can be installed and lauched as a Python Dash application on your local machine. This procedure requires prior installation of [Anaconda](https://www.anaconda.com/) or [miniconda3](https://docs.conda.io/en/latest/miniconda.html).

### Step 1. Create Environment:
Clone this repository, navigate to the goloco main directory, and run the following commands to create and activate a new python virtual environment:

```bash
conda create --name goloco -c conda-forge python=3.11
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
goloco uses a modification of the chronos software by the Broad Institute to generate psuedo-chronos gene effect scores from smaller scale CRISPR knock out experiments (i.e, lossy subset experiments). Visit the [chronos](https://github.com/broadinstitute/chronos) documentation for any specific installation instructions than may be relevant to your machine.

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
Prior to lauching this application, please see the data section for instructions on how to download the required models to your local hard drive and to modify application configurations and source code to access these files.

This app uses three main services to operate--Redis, Celery, and Python Dash. To lauch these services open three terminals in the main goloco directory, activate the virtual environment in each, and run the following commands: 

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
goloco should now be operational and serviced over port 8080 on your local machine. Open any web browser and type '127.0.0.1:8080' or 'localhost:8080' to access the application. Note that the port can be changed in the app.py source code if necessary or desired. The application can be stopped by closing processes on each of the service terminals and restarted following Step 4.

## Data:



## References

[1] 
