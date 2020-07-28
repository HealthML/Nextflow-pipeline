# NGS580 container for making reports with R 3.4.3
FROM ubuntu:18.04

MAINTAINER Alva James

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
    apt-get install -y wget \
    bzip2 \
    tabix \
    unzip  \
    curl  \
    gcc  \
    g++  \
    build-essential  \
    python-dev \
    libevent-dev  \
    git

#RUN mkdir testdir/
#RUN git clone https://github.com/HealthML/seak.git testdir


# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
    bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
    rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ VEP ~~~~~ #
RUN conda install -y -c bioconda \
    ensembl-vep=99.0 \
    gatk4=4.0.5.1 \
    numpy \
    pandas \
    cython

RUN git clone https://github.com/HealthML/seak.git
RUN pip install -e ./seak

COPY run_test_proteinlof.py ./run_test_proteinlof.py
COPY filter_VEP.py ./filter_VEP.py
RUN find . -name "*seak*" -print
RUN find . -name "*VEP.py" -print
RUN realpath filter_VEP.py
RUN find $PWD -type f | grep "filter_VEP.py"



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

#COPY /home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak seak/
# The <src> path must be inside the context of the build; you cannot ADD ../something /something, because the first step of a docker build is to send the context directory (and subdirectories) to the docker daemon.ss
# Multiple <src> resources may be specified but the paths of files and directories will be interpreted as relative to the source of the context of the build. (relative to the source directory that is being built)
