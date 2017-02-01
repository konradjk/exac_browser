#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Bring down previous replication controller
# kubectl delete service exac-serve
kubectl delete service mongo
kubectl delete rc mongo-controller
kubectl delete pod exac-load
# kubectl delete rc exac-serve


# Delete the cluster
# gcloud container clusters delete exac-loading-cluster --zone us-east1-d