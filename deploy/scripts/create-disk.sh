#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Create the disc
gcloud compute disks create --size=200GB --zone=us-east1-d mongo-disk

# Delete
# gcloud compute disks delete mongo-disk --zone=us-east1-d 