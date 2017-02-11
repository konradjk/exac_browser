#!/bin/bash

# halt on any error
set -e

gcloud container clusters create gnomad-serving-cluster \
--machine-type n1-standard-4 \
--zone us-east1-d \
--num-nodes 1 \
--project exac-gnomad
