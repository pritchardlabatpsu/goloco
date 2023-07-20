FROM python:3.11

RUN set -ex \
    && pip install --upgrade pip \
    && apt-get update

ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

RUN pip install --no-cache-dir -r requirements.txt --use-pep517
RUN pip install scikit-learn==0.22.1
RUN git clone https://github.com/broadinstitute/chronos.git
RUN python /app/chronos/setup.py install
RUN unzip $(find /usr/local/lib/python3.11/site-packages -name '*.egg') -d /usr/local/lib/python3.11/site-packages

#ENV PORT=8080

#EXPOSE $PORT

#CMD gunicorn app:server -b 0.0.0.0:${PORT}