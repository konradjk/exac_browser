#!/bin/bash

# halt on any error
set -e

# Set project and kubernetes cluster
gcloud config set project exac-gnomad
kubectl config use-context gke_exac-gnomad_us-east1-d_exac-serving-cluster

# Start mongo -- takes 20 seconds or so
kubectl create -f deploy/config/mongo-service.yaml
kubectl create -f deploy/config/mongo-controller.yaml
sleep 30

# Start the server and expose to the internet w/ autoscaling & load balancing
kubectl create -f deploy/config/exac-serve-rc.json
kubectl expose rc exac-serve --type="LoadBalancer" --load-balancer-ip=35.185.33.81
kubectl autoscale rc exac-serve --min=1 --max=2 --cpu-percent=80
