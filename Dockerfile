FROM ubuntu:22.04

#Installing python

RUN apt update
RUN apt -y upgrade

RUN apt install -y python3-pip

RUN pip install dna_features_viewer==3.1.1

RUN pip install svg_stack==0.1.0

RUN pip install biopython==1.79

RUN pip install matplotlib==3.3.4

#To manage the time zone required when installing R

ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

#Installing R

RUN apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt install -y r-base

RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev

RUN Rscript -e 'install.packages("ggplot2")'

RUN Rscript -e 'install.packages("tidyverse")'

#prodigal & phanotate

RUN pip install phanotate

RUN apt-get install prodigal=1:2.6.3-5

#samtools 

RUN apt install wget

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2

RUN tar xvjf samtools-1.10.tar.bz2

RUN cd samtools-1.10
RUN make .

ENV PATH=$PATH:/samtools-1.10

RUN cd /

#bcftools

RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2

RUN tar xvjf bcftools-1.14.tar.bz2

RUN cd bcftools-1.14
RUN make .

ENV PATH=$PATH:/bcftools-1.14

RUN cd /

#bwa
RUN apt-get install bwa=0.7.17-6

#codingDiv repository

RUN apt install -y git

RUN git clone https://github.com/ericolo/potential-garbanzo

RUN chmod +x /potential-garbanzo/garbanzo/*

ENV PATH=$PATH:/potential-garbanzo/garbanzo

WORKDIR /data

