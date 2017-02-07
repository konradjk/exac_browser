#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the disks
gcloud compute disks create --size=200GB --zone=us-east1-d mongo-disk
gcloud compute disks create --size=200GB --zone=us-east1-d exac-logs

# Delete disks
# gcloud compute disks delete mongo-disk --zone=us-east1-d
# gcloud compute disks delete --size=200GB --zone=us-east1-d exac-logs