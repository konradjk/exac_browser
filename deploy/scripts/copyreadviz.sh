#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the disk, wait 5 minutes
gcloud compute disks create --size=1000GB --zone=us-east1-d exac-readviz-exons-gpd

"$(dirname "$0")"/start-load-cluster.sh
kubectl config use-context gke_exac-gnomad_us-east1-d_gnomad-loading-cluster

# docker build -f deploy/dockerfiles/copyreadviz.dockerfile -t gcr.io/exac-gnomad/copyreadviz .
# gcloud docker push gcr.io/exac-gnomad/copyreadviz

sleep 60

# start pod
kubectl create -f deploy/config/copyreadviz-pod.json

# takedown pod on finish (about an hour?)
# kubectl delete pod exac-copyreadviz

# Delete
# gcloud compute disks delete exacv1-readviz --zone=us-east1-d

# takedown cluster
# "$(dirname "$0")"/takedown-load.sh
