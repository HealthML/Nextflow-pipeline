# NGS580 container for making reports with R 3.4.3
FROM ubuntu:18.04

MAINTAINER Alva James

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
    apt-get install -y wget \
    bzip2 \
    tabix \
    unzip  \
    curl

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
    bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
    rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ VEP ~~~~~ #
RUN conda install -y -c bioconda \
    ensembl-vep=99.0 \
    gatk4=4.0.5.1

ENV PLINK_VERSION v2.00a2.3LM
ENV PLINK_ZIP plink_linux_x86_64.zip
ENV PLINK_DIR /opt/plink2

# ~~~~~~ PLINK2 ~~~~#
RUN curl http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip -o $PLINK_ZIP \
    && mkdir $PLINK_DIR \
    && unzip $PLINK_ZIP -d $PLINK_DIR \
    && rm $PLINK_ZIP

RUN cd /usr/local/bin \
    && ln -s $PLINK_DIR/plink2 plink2


# We will define plink2 as the application to run when this container starts
ENTRYPOINT ["plink2"]

