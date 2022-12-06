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

#Installing specific versions to avoid update problems
RUN apt install wget

#tidyverse
RUN wget http://cran.r-project.org/src/contrib/Archive/tidyverse/tidyverse_1.3.0.tar.gz

#To get all dependencies
RUN Rscript -e 'install.packages("tidyverse")'

#Then the right version
RUN R CMD INSTALL tidyverse_1.3.0.tar.gz

#ggplot2
RUN wget http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.5.tar.gz

RUN Rscript -e 'install.packages("ggplot2")'

RUN R CMD INSTALL ggplot2_3.3.5.tar.gz

#prodigal & phanotate & getorf

RUN pip install phanotate

RUN apt-get install prodigal=1:2.6.3-5

RUN apt-get install -y emboss=6.6.0+dfsg-11ubuntu1

#samtools 

RUN bash

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2

RUN tar xvjf samtools-1.10.tar.bz2

RUN make -C samtools-1.10

ENV PATH=$PATH:/samtools-1.10

RUN echo $PATH

#bcftools

RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2

RUN tar xvjf bcftools-1.14.tar.bz2

RUN make -C bcftools-1.14

ENV PATH=$PATH:/bcftools-1.14

RUN echo $PATH

#bwa
RUN apt-get install bwa=0.7.17-6

#svglite to save snp plot as svg

RUN apt install -y libfontconfig1-dev

RUN Rscript -e 'install.packages("svglite")'

RUN pip install lxml

#codingDiv repository

RUN apt install -y git

#To not cach the git clone command as repository can change
#This will change if the github repo is modified
ADD https://api.github.com/repos/ericolo/codingDiv/git/refs/heads/main version.json

RUN git clone https://github.com/ericolo/codingDiv

RUN chmod +x /codingDiv/scripts/*

ENV PATH=$PATH:/codingDiv/scripts

RUN echo $PATH

WORKDIR /data

#To not type codingDiv.sh 
CMD ["bash", "codingDiv.sh"]

#It works like this, files will be written in the given dir
#docker run -v /Users/ONE/Downloads/codingdiv:/data codingdiv codingDiv.sh tylcv.fna blast_hits_90.fna 90 1 2 1 3 N

#docker run -v /Users/ONE/Downloads/codingdiv:/data codingdiv tylcv.fna blast_hits_90.fna 90 1 2 1 3 N

#docker build --tag codingdiv .
