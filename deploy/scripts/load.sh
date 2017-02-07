#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the replication controller
kubectl create -f deploy/config/mongo-service.yaml
kubectl create -f deploy/config/mongo-controller.yaml

# Wait for mongo to initialize
sleep 120

# load data
kubectl create -f deploy/config/exac-load-pod.json

# copy readviz files
kubectl create -f deploy/config/readviz-pod.json
kubectl create -f deploy/config/readviz-pod.json
