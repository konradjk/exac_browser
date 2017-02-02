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
sleep 30
kubectl create -f deploy/config/exac-serve-rc.json

# Expose replication controller to the internet
kubectl expose rc exac-serve --type="LoadBalancer"
kubectl autoscale rc exac-serve --min=5 --max=20 --cpu-percent=30