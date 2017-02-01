#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
gcloud container clusters create exac-serving-cluster --machine-type n1-standard-4 --zone us-east1-d --num-nodes 2
