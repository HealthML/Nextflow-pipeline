# NGS580 container for making reports with R 3.4.3
FROM ubuntu:18.04

MAINTAINER Alva James

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ VEP ~~~~~ #
RUN conda install -y -c bioconda \
ensembl-vep=99.0 \
gatk4=4.0.5.1

# ~~~~~ OTHER PYTHON PACKAGES ~~~~~ #
RUN conda install --yes \
    numpy \
    pandas \
    cython

#COPY /home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak seak/
# The <src> path must be inside the context of the build; you cannot ADD ../something /something, because the first step of a docker build is to send the context directory (and subdirectories) to the docker daemon.
