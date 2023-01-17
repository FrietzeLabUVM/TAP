FROM ubuntu:focal
RUN apt-get update
#needed to install git without prompts
RUN DEBIAN_FRONTEND=noninteractive TZ=America/New_York apt-get -y install tzdata
RUN apt-get install git-all -y
## samtools deps
RUN apt-get install libcurl4-gnutls-dev -y
RUN apt-get install libssl-dev -y
RUN apt-get install libncurses-dev -y
RUN apt-get install zlib1g-dev -y
RUN apt-get install libbz2-dev -y
RUN apt-get install liblzma-dev -y
# Install base utilities
RUN apt-get install -y build-essential  && \
    apt-get install -y wget 
## samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
  tar -xf samtools-1.16.1.tar.bz2 && \
  cd samtools-1.16.1 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install
## python
RUN apt-get install software-properties-common -y && \
  apt-add-repository ppa:deadsnakes/ppa -y && \
  apt-get install python3.8 -y && \
  apt-get install python3-pip -y

## bedtools
# RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
#   tar -xzf bedtools-2.30.0.tar.gz && \
#   cd bedtools2/ && \
#   make && \
#   mv bin/* /usr/local/bin/
#RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
#  mv bedtools.static.binary /usr/bin/bedtools && \
#  chmod a+x /usr/bin/bedtools
RUN apt-get install bedtools

## bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
  tar -xjf bcftools-1.10.2.tar.bz2 && \
  cd bcftools-1.10.2 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install

## STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
  tar -xzf 2.7.10b.tar.gz && \
  cd STAR-2.7.10b/source && \
  make STAR && \
  ln -s /STAR-2.7.10b/source/STAR /usr/local/bin
## MACS2
#RUN conda install -c "bioconda/label/cf201901" macs2
RUN pip install macs3 && \
  ln -s /usr/local/bin/macs3 /usr/local/bin/macs2

