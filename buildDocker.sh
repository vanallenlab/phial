#!/bin/bash

tag=1.0.3

# Replace ./ with /./ for Docker usage on FireCloud
sed -i '' "s#Rscript PHIAL_v1.0.R#Rscript /PHIAL_v1.0.R#" PHIAL_wrapper.R
sed -i '' "s#Nozzle_template.R#/Nozzle_template.R#" PHIAL_v1.0.R
sed -i '' "s#empty_inputs/#/empty_inputs/#" PHIAL_v1.0.R

docker build -t vanallenlab/phial:${tag} .
docker push vanallenlab/phial:${tag}

# Replace /./ with ./ to return to local usage
sed -i '' "s#Rscript /PHIAL_v1.0.R#Rscript PHIAL_v1.0.R#" PHIAL_wrapper.R
sed -i '' "s#/Nozzle_template.R#Nozzle_template.R#" PHIAL_v1.0.R
sed -i '' "s#/empty_inputs/#empty_inputs/#" PHIAL_v1.0.R
