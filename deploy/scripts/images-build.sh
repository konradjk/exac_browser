#!/bin/bash

# halt on any error
set -e

# Build docker images
docker build -f deploy/dockerfiles/exacbase.dockerfile -t gcr.io/exac-gnomad/exacbase .
docker build -f deploy/dockerfiles/exacserve.dockerfile -t gcr.io/exac-gnomad/exacserve .
docker build -f deploy/dockerfiles/exacload.dockerfile -t gcr.io/exac-gnomad/exacload .
