#!/bin/bash

# halt on any error
set -e

# Set project
gcloud config set project exac-gnomad

# Push docker images
# gcloud docker push gcr.io/exac-gnomad/exacbase
gcloud docker push gcr.io/exac-gnomad/exaccode
gcloud docker push gcr.io/exac-gnomad/exacload
# gcloud docker push gcr.io/exac-gnomad/exacserve
