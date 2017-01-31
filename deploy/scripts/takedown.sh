#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Bring down previous replication controller
kubectl delete service exac-controller
kubectl delete service mongo
kubectl delete rc mongo-controller
kubectl delete rc exac-controller

# Delete the cluster
gcloud container clusters delete exac-container-cluster --zone us-east1-d