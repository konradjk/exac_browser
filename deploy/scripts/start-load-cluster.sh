#!/bin/bash

# halt on any error
set -e

gcloud container clusters create exac-loading-cluster \
--machine-type n1-highmem-8 \
--zone us-east1-d \
--num-nodes 1 \
--project exac-gnomad
