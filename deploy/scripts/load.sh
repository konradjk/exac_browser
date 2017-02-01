#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
# gcloud container clusters create exac-container-cluster --machine-type f1-micro --zone us-east1-d

# Create the replication controller
kubectl create -f deploy/config/mongo-service.yaml
kubectl create -f deploy/config/mongo-controller.yaml

# Wait for mongo to initialize
sleep 120

# load data
kubectl create -f deploy/config/exac-load-pod.json
