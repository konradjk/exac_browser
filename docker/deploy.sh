#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
gcloud container clusters create exac-container-cluster --machine-type f1-micro --zone us-east1-d

# Build docker images
# docker build -f docker/exacserve.dockerfile -t gcr.io/exac-gnomad/exacserve:v1 .
docker build -f docker/exacbase.dockerfile -t macarthurlab/exacbase .
docker build -f docker/exacserve.dockerfile -t macarthurlab/exacserve .
docker build -f docker/exacload.dockerfile -t macarthurlab/exacload .

# Push docker images
# gcloud docker push gcr.io/exac-gnomad/exacserve
docker push macarthurlab/exacbase
docker push macarthurlab/exacserve
docker push macarthurlab/exacload

# Bring down previous replication controller
kubectl delete service exac-controller
kubectl delete service mongo
kubectl delete rc mongo-controller
kubectl delete rc exac-controller

# Create the replication controller
kubectl create -f mongo-service.yaml
kubectl create -f mongo-controller.yaml
kubectl create -f kubernetes.json

# Expose replication controller to the internet
kubectl expose rc exac-controller --type="LoadBalancer"
