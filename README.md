# docker_4dn
Repository of Docker image for 4dn pipeline. Docker image is at: https://hub.docker.com/repository/docker/rnakato/4dn/general

## Run

For Docker:

    # pull docker image
    docker pull rnakato/4dn

    # container login
    docker run [--gpus all] --rm -it rnakato/4dn /bin/bash
    # jupyter notebook
    docker run [--gpus all] --rm -p 8888:8888 -v (your directory):/opt/work rnakato/4dn <command>

For Singularity:

    # build image
    singularity build -F rnakato_4dn.sif docker://rnakato/4dn 
    # jupyter notebook
    singularity exec [--nv] rnakato_4dn.sif <command>

## Build image from Dockerfile
First clone and move to the repository

    git clone https://github.com/rnakato/docker_4dn.git
    cd docker_4dn

Then type:

    docker build -t <account>/4dn .

## Contact

Ryuichiro Nakato: rnakato AT iqb.u-tokyo.ac.jp
