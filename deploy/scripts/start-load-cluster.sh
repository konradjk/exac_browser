#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
# gcloud container clusters create exac-container-cluster --machine-type f1-micro --zone us-east1-d
gcloud container clusters create exac-loading-cluster --machine-type n1-standard-8 --zone us-east1-d --num-nodes 1
# gcloud container clusters create exac-loading-cluster --machine-type n1-standard-32 --zone us-east1-d --num-nodes 1
# gcloud container clusters create exac-loading-cluster --machine-type n1-highmem-32 --zone us-east1-d --num-nodes 1
