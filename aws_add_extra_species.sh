#!/bin/bash
## Set AWS credentials
export AWS_ACCESS_KEY_ID=AKIAWHC4XC6G3CEJ5ZXJ
export AWS_SECRET_ACCESS_KEY=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
export DEFAULT_REGION_NAME=us-east-2
export DEFAULT_OUTPUT_FORMAT=json

#download and unzip orthofinder
wget https://github.com/davidemms/OrthoFinder/releases/download/2.5.1/OrthoFinder.tar.gz
tar xzvf OrthoFinder.tar.gz

aws s3 cp s3://mukherjee-lab/dash_app_pseudomonas_1/2023-07-27_OF_Results.tar.gz ./
aws s3 cp s3://mukherjee-lab/dash_app_pseudomonas_2/New_proteomes/ ./New_proteomes --recursive

#download yum and decompress results folder
sudo yum install pigz -y
pigz -dc 2023-07-27_OF_Results.tar.gz | tar xv -

./OrthoFinder/orthofinder  -b ./2023-07-27_OF_Results/Results_Jul27 -f ./New_proteomes


#move species tree to s3
aws s3 cp ./2023-07-27_OF_Results/Results_Aug04/Species_Tree/SpeciesTree_rooted_node_labels.txt s3://mukherjee-lab/dash_app_pseudomonas_2/


time tar -I pigz -cf Results_Aug04.tar.gz 2023-07-27_OF_Results/Results_Aug04

aws s3 cp Results_Aug04.tar.gz s3://mukherjee-lab/dash_app_pseudomonas_2/


# if new tree is used:
./OrthoFinder/orthofinder  -ft 2023-07-27_OF_Results/Results_Aug04/ -s ps_reroot.txt