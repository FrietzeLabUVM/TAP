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

## salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz && \
  tar -xf salmon-1.9.0_linux_x86_64.tar.gz && \
  cp salmon-1.9.0_linux_x86_64/bin/salmon /usr/local/bin/ && \
  cp -n salmon-1.9.0_linux_x86_64/lib/* /usr/lib/x86_64-linux-gnu

## suppa2
RUN pip install SUPPA==2.3 && \
  git clone https://github.com/comprna/SUPPA && \
  sed -i '1s/^/\#\!\/usr\/bin\/python3 \n/' SUPPA/suppa.py && \
  chmod a+x SUPPA/*py && \
  chmod a+x SUPPA/*pl && \
  cp -n SUPPA/*py /usr/local/bin/  && \
  cp -n SUPPA/*pl /usr/local/bin/ 
  
#RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedSort  && \
#  chmod a+x bedSort   && \
#  mv bedSort /usr/local/bin/

## bedGraphToBigWig
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
  chmod a+x bedGraphToBigWig  && \
  mv bedGraphToBigWig /usr/local/bin/

## All UCSC, too big 3.8G
# RUN apt-get install rsync -y && \
#   mkdir ucsc && \
#   cd ucsc && \
#   rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./ && \
#   chmod a+x * && \
#   mv * /usr/local/bin/

## exactSNP
RUN wget -N https://sourceforge.net/projects/subread/files/subread-2.0.3/subread-2.0.3-source.tar.gz && \
  tar -xf subread-2.0.3-source.tar.gz  && \
  cd subread-2.0.3-source/src && \
  make -f Makefile.Linux && \
  cd ../bin/ && \
  mv * /usr/local/bin/

## pigz
RUN wget http://zlib.net/pigz/pigz-2.7.tar.gz && \
  tar -xf pigz-2.7.tar.gz && \
  cd pigz-2.7 && \
  make && \
  mv pigz /usr/local/bin/ && \
  mv unpigz /usr/local/bin/ && \
  cd .. && \
  rm -r pigz-2.7

## permissions and cleanup
RUN chown -R root:root salmon-1.9.0_linux_x86_64 && \
  chown -R root:root subread-2.0.3-source && \
  rm *tar*

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-ubuntu64.tar.gz && \
  tar -xf sratoolkit.3.0.5-ubuntu64.tar.gz && \
  mv sratoolkit.3.0.5-ubuntu64/bin/* /usr/local/bin/ && \
  rm -r sratoolkit.3.0.5-ubuntu64
