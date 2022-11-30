#!/bin/bash

# Build the docker image with the tag "py37_granuloma_scrna" and the same user name and group permission as the host user
if [[ $OSTYPE == 'darwin'* ]]; then
    echo 'detected MacOS system'
    docker build Docker -f Docker/Dockerfile.macos -t py37_granuloma_scrnaseq --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) 
elif [[ $OSTYPE == 'linux-gnu' ]]; then
    echo 'detected Linux system'
    docker build Docker -f Docker/Dockerfile.linux -t py37_granuloma_scrnaseq --build-arg USER_ID=$(id -u) --build-arg GROUP_ID=$(id -g) 
else
    echo 'not a supported OS. Please use MacOS or Linux to build this image.'
fi
