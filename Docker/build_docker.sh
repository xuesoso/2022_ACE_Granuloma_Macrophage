#!/bin/bash

# Build the docker image with the tag "py37_granuloma_scrna" and the same user name and group permission as the host user
docker build Docker -f Docker/Dockerfile -t py37_granuloma_scrnaseq --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) 
