#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
gcloud container clusters create exac-container-cluster --machine-type f1-micro --zone us-east1-d

# Create the replication controller
kubectl create -f ../config/mongo-service.yaml
kubectl create -f ../config/mongo-controller.yaml
kubectl create -f ../config/kubernetes.json

# Expose replication controller to the internet
kubectl expose rc exac-controller --type="LoadBalancer"
