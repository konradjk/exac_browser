FROM python:2.7.12

MAINTAINER MacArthur Lab

ENV EXAC_DATA=/var/exac_data/ READ_VIZ=/mongo/readviz \
  CLOUD_SDK_REPO=cloud-sdk-jessie GCSFUSE_REPO=gcsfuse-jessie

COPY . /var/www
WORKDIR /var/www

RUN apt-get update && apt-get install -y apt-transport-https

# Install gcloud
RUN echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
  apt-get update && apt-get install -y google-cloud-sdk

# Install gcfuse
RUN echo "deb http://packages.cloud.google.com/apt $GCSFUSE_REPO main" | tee /etc/apt/sources.list.d/gcsfuse.list && \
  curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
  apt-get update && apt-get -y install gcsfuse

# Authorize container with service account
RUN gcloud auth activate-service-account \
  --key-file=deploy/keys/exac-gnomad-30ea80400948.json

RUN pip install -r requirements.txt

