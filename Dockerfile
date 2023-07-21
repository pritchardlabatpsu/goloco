FROM python:3.8

RUN set -ex \
    && pip install --upgrade pip \
    && apt-get update

ENV APP_HOME /app
WORKDIR $APP_HOME
COPY . ./

RUN pip install --no-cache-dir -r requirements.txt --use-pep517
RUN pip install scikit-learn==0.22.1
RUN git clone https://github.com/broadinstitute/chronos.git
RUN cd /app/chronos \
    && pip install .

#ENV PORT=8080

#EXPOSE $PORT

#CMD gunicorn app:server -b 0.0.0.0:${PORT}