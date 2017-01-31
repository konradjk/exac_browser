#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the cluster
gcloud container clusters create exac-container-cluster --machine-type f1-micro --zone us-east1-d
