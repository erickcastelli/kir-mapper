FROM ubuntu:22.04

WORKDIR /home

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV DEBIAN_FRONTEND noninteractive
ENV TZ UTC


RUN apt-get update && apt-get install -y build-essential cmake zlib1g-dev libncurses-dev git libboost-all-dev wget libcurl4-openssl-dev gzip zip libssl-dev nano

RUN apt-get install -y python3-dev python3-pip

RUN pip3 install --user whatshap==2.2
RUN export PATH=$HOME/.local/bin:$PATH

RUN apt-get install -y --no-install-recommends \
                 littler \
 		 r-base \
 		 r-base-dev \
 		 r-recommended \
         r-cran-docopt
		 
RUN apt-get install -y --no-install-recommends \
                 bwa=0.7.17-6


RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2
RUN tar -xvf samtools-1.19.2.tar.bz2
WORKDIR /home/samtools-1.19.2
RUN make
RUN make install

WORKDIR /home
RUN wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
RUN tar -xvf bcftools-1.19.tar.bz2
WORKDIR /home/bcftools-1.19
RUN make
RUN make install

WORKDIR /home
RUN mkdir freebayes
WORKDIR /home/freebayes
RUN wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz
RUN gunzip freebayes-1.3.6-linux-amd64-static.gz
RUN chmod 777 freebayes-1.3.6-linux-amd64-static
RUN mv freebayes-1.3.6-linux-amd64-static freebayes
RUN cp freebayes /usr/bin/


WORKDIR /home
RUN wget https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip
RUN unzip STAR_2.7.11b.zip
WORKDIR /home/STAR_2.7.11b
RUN cp Linux_x86_64_static/STAR /usr/bin


WORKDIR /home
RUN mkdir shapeit4
WORKDIR /home/shapeit4
RUN wget https://github.com/odelaneau/shapeit4/archive/refs/tags/v4.2.2.tar.gz
RUN tar -xvf v4.2.2.tar.gz
WORKDIR /home/shapeit4/shapeit4-4.2.2

WORKDIR /home
RUN mkdir Tools
WORKDIR /home/Tools
RUN wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2
RUN tar -xvf htslib-1.11.tar.bz2
WORKDIR /home/Tools/htslib-1.11
RUN make

WORKDIR /home/shapeit4/shapeit4-4.2.2
RUN make HOME=/home/
WORKDIR /home/shapeit4/shapeit4-4.2.2/bin
RUN mv shapeit4.2 shapeit4
RUN cp shapeit4 /usr/bin


WORKDIR /home

RUN apt-get install -y openjdk-18-jdk
	
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('plotly')"
RUN Rscript -e "install.packages('htmlwidgets')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('pacman')"
RUN Rscript -e "install.packages('shiny')"
RUN Rscript -e "install.packages('shinyjs')"
RUN Rscript -e "install.packages('purrr')"
RUN Rscript -e "install.packages('readr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('tidyr')"
RUN Rscript -e "install.packages('DT')"
RUN Rscript -e "install.packages('forcats')"
RUN Rscript -e "install.packages('shinyalert')"


RUN wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar


WORKDIR /home
RUN git clone https://github.com/erickcastelli/kir-mapper
WORKDIR /home/kir-mapper
RUN mkdir /home/kir-mapper/build
WORKDIR /home/kir-mapper/build
RUN cmake ../src/
RUN make
RUN cp kir-mapper /usr/bin


RUN cp ~/.local/bin/whatshap /usr/bin/



WORKDIR /home
RUN echo "bwa=/usr/bin/bwa" >> ~/.kir-mapper
RUN echo "samtools=/usr/local/bin/samtools" >> ~/.kir-mapper
RUN echo "bcftools=/usr/local/bin/bcftools" >> ~/.kir-mapper
RUN echo "freebayes=/usr/bin/freebayes" >> ~/.kir-mapper
RUN echo "shapeit4=/usr/bin/shapeit4" >> ~/.kir-mapper
RUN echo "star=/usr/bin/STAR" >> ~/.kir-mapper
RUN echo "picard=/home/picard.jar" >> ~/.kir-mapper
RUN echo "whatshap=/usr/bin/whatshap" >> ~/.kir-mapper
RUN echo "db=/home/kir-mapper/kir-mapper_db_latest" >> ~/.kir-mapper

WORKDIR /home/kir-mapper
RUN wget --no-check-certificate https://www.castelli-lab.net/support/kir-mapper_db_latest.zip
RUN unzip kir-mapper_db_latest.zip
RUN rm kir-mapper_db_latest.zip

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y locales

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG en_US.UTF-8


RUN apt-get install -y tabix

RUN apt-get install -y parallel

RUN apt-get install -y htop



WORKDIR /home

CMD ["bash"]


