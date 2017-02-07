#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the replication controller
kubectl config use-context gke_exac-gnomad_us-east1-d_exac-serving-cluster
kubectl create -f deploy/config/mongo-service.yaml
kubectl create -f deploy/config/mongo-controller.yaml
sleep 30
kubectl create -f deploy/config/exac-serve-rc.json

# Expose replication controller to the internet
kubectl expose rc exac-serve --type="LoadBalancer" --load-balancer-ip=35.185.33.81
kubectl autoscale rc exac-serve --min=1 --max=2 --cpu-percent=80
