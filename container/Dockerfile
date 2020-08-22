# NGS580 container for making reports with R 3.4.3
FROM ubuntu:18.04

MAINTAINER Alva James

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
    apt-get install -y wget \
    bzip2 \
    tabix \
    unzip  \
    curl \
    gcc  \
    g++  \
    build-essential  \
    python-dev \
    libevent-dev \
    git  
    

# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
    bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
    rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"

# ~~~~~ VEP & SEAK-libraries~~~~~ #
RUN conda install -y -c bioconda \
    ensembl-vep=99.0 \
    gatk4=4.0.5.1

ENV PLINK_VERSION v2.00a2.3LM
ENV PLINK_ZIP plink_linux_x86_64.zip
ENV PLINK_DIR /opt/plink2

# ~~~~~~ PLINK2 ~~~~#
RUN curl http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip -o $PLINK_ZIP \
    && mkdir $PLINK_DIR \
    gatk4=4.0.5.1 \
    numpy \
    pandas \
    cython

# ~~~~~~~~~~~~~~~~~ SEAK ~~~~~~~~~~ #
#RUN apt-get install git
RUN git clone https://github.com/HealthML/seak.git
RUN pip install -e ./seak

COPY run_test_proteinlof.py ./run_test_proteinlof.py
COPY filter_VEP.py ./filter_VEP.py
RUN find . -name "*seak*" -print
RUN find . -name "*VEP.py" -print
RUN realpath filter_VEP.py
RUN find $PWD -type f | grep "filter_VEP.py"

#RUN pip install graphviz
#FROM alpine:3.7
#RUN apk update && apk add bash && apk add coreutils && apk add curl && apk add graphviz

#ENTRYPOINT "/usr/bin/curl"

# ~~~~~~ PLINK2 ~~~~#
ENV PLINK_VERSION v2.00a2.3LM
ENV PLINK_ZIP plink_linux_x86_64.zip
ENV PLINK_DIR /opt/plink2

RUN curl http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip -o $PLINK_ZIP \
    && mkdir -p $PLINK_DIR \
    && unzip $PLINK_ZIP -d $PLINK_DIR \
    && rm $PLINK_ZIP

RUN cd /usr/local/bin \
    && ln -s $PLINK_DIR/plink2 plink2


=======
# We will define plink2 as the application to run when this container starts
ENTRYPOINT ["plink2"]

