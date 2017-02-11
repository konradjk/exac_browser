#!/bin/bash

# halt on any error
set -e

gcloud container clusters create gnomad-loading-cluster \
--machine-type n1-highmem-32 \
--zone us-east1-d \
--num-nodes 1 \
--project exac-gnomad
