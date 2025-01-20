FROM python:3.9

ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

RUN apt-get update && apt-get install -y libhdf5-dev
RUN pip install --upgrade pip
RUN apt-get install -y pkg-config

RUN pip install python-dotenv
RUN pip install --no-binary h5py h5py
RUN pip install --no-cache-dir -r requirements.txt --use-pep517
RUN cd /app/chronos && pip install .
RUN pip install scikit-learn==0.22.1