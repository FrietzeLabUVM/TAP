FROM ubuntu:focal-20230605
RUN apt-get update

## Install base utilities and deps
## needed to install git without prompts
RUN DEBIAN_FRONTEND=noninteractive TZ=America/New_York apt-get -y install tzdata && \
  apt-get install git-all -y && \
  apt-get install libcurl4-gnutls-dev -y && \
  apt-get install libssl-dev -y && \
  apt-get install libncurses-dev -y && \
  apt-get install zlib1g-dev -y && \
  apt-get install libbz2-dev -y && \
  apt-get install liblzma-dev -y && \
  apt-get install -y build-essential  && \
  apt-get install -y wget 


#apt-get install software-properties-common -y && \
# ## python
RUN apt-get install software-properties-common -y && \
  apt-add-repository ppa:deadsnakes/ppa -y && \
  apt-get install python3.11 python3.11-dev -y && \
  ln -sf /usr/bin/python3.11 /usr/bin/python3 && \
  wget https://bootstrap.pypa.io/get-pip.py && \
  python3 get-pip.py
  
## MACS2
RUN pip install macs3==3.0.0b3 && \
  ln -s /usr/local/bin/macs3 /usr/local/bin/macs2

## samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
  tar -xf samtools-1.16.1.tar.bz2 && \
  cd samtools-1.16.1 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && \
  rm -r samtools-1.16.1 && \
  rm samtools-1.16.1.tar.bz2 

## bedtools
RUN apt-get install bedtools

## bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 && \
  tar -xjf bcftools-1.10.2.tar.bz2 && \
  cd bcftools-1.10.2 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && \
  rm -r bcftools-1.10.2 && \
  rm bcftools-1.10.2.tar.bz2

## STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
  tar -xzf 2.7.10b.tar.gz && \
  cd STAR-2.7.10b/source && \
  make STAR && \
  cd ../.. && \
  mv STAR-2.7.10b/source/STAR /usr/local/bin && \
  rm -r STAR-2.7.10b && \
  rm 2.7.10b.tar.gz



## salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
  tar -xf salmon-1.9.0_linux_x86_64.tar.gz && \
  cp salmon-1.9.0_linux_x86_64/bin/salmon /usr/local/bin/ && \
  cp -n salmon-1.9.0_linux_x86_64/lib/* /usr/lib/x86_64-linux-gnu && \
  rm -r salmon-1.9.0_linux_x86_64 && \
  rm salmon-1.9.0_linux_x86_64.tar.gz

## suppa2
RUN pip install SUPPA==2.3 && \
  git clone https://github.com/comprna/SUPPA && \
  sed -i '1s/^/\#\!\/usr\/bin\/python3 \n/' SUPPA/suppa.py && \
  chmod a+x SUPPA/*py && \
  chmod a+x SUPPA/*pl && \
  cp -n SUPPA/*py /usr/local/bin/ && \
  cp -n SUPPA/*pl /usr/local/bin/ && \
  rm -r SUPPA

## bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
  chmod a+x bedGraphToBigWig  && \
  mv bedGraphToBigWig /usr/local/bin/

## exactSNP
RUN wget -N https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-source.tar.gz && \
  tar -xf subread-2.0.3-source.tar.gz --no-same-owner && \
  cd subread-2.0.3-source/src && \
  make -f Makefile.Linux && \
  cd ../bin/ && \
  mv * /usr/local/bin/ && \
  cd ../.. && \
  rm -r subread-2.0.3-source && \
  rm subread-2.0.3-source.tar.gz

## pigz
RUN wget http://zlib.net/pigz/pigz-2.8.tar.gz && \
  tar -xf pigz-2.8.tar.gz && \
  cd pigz-2.8 && \
  make && \
  mv pigz /usr/local/bin/ && \
  mv unpigz /usr/local/bin/ && \
  cd .. && \
  rm -r pigz-2.8 && \
  rm pigz-2.8.tar.gz

## SRA toolkit, fasterqdump
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-ubuntu64.tar.gz && \
  tar -xf sratoolkit.3.0.5-ubuntu64.tar.gz && \
  mv sratoolkit.3.0.5-ubuntu64/bin/* /usr/local/bin/ && \
  rm -r sratoolkit.3.0.5-ubuntu64 && \
  rm sratoolkit.3.0.5-ubuntu64.tar.gz

## bedops, provides RAM limited alternative to bedSort
RUN mkdir bedops  && \
  cd bedops  && \
  wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 && \
  tar -xf bedops_linux_x86_64-v2.4.41.tar.bz2 && \
  cp bin/sort-bed /usr/local/bin/ && \
  cp bin/sort-bed-typical /usr/local/bin/ && \
  cd .. && \
  rm -r bedops

## gffread, used to setup references
RUN git clone https://github.com/gpertea/gffread  && \
  cd gffread  && \
  git reset --hard 558c97f && \
  make release && \
  mv gffread /usr/local/bin/ && \
  cd .. && \
  rm -r gffread