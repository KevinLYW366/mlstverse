FROM ubuntu:22.04

MAINTAINER matsumoto <matsumoto@gen-info.osaka-u.ac.jp>

ENV DEBIAN_FRONTEND=noninteractive
ENV VERSION 0.1.4-46485e3

RUN apt-get update \
 && apt-get -y install --no-install-recommends \
    automake \
    autoconf \
    build-essential \
    ca-certificates \
    curl \
    git \
    libcairo2-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libpcre2-dev \
    libpng-dev \
    libssl-dev \
    libtiff5-dev \
    libxml2-dev \
    pkg-config \
    r-base \
    r-base-dev \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/src

RUN git clone https://github.com/lh3/minimap2 \
 && cd minimap2 \
 && make \
 && cp minimap2 /usr/bin/

RUN git clone https://github.com/lh3/bwa \
 && cd bwa \
 && make \
 && cp bwa /usr/bin/

RUN git clone --recursive https://github.com/samtools/htslib \
 && cd htslib \
 && autoreconf -fi \
 && ./configure \
 && make \
 && make install

RUN git clone https://github.com/samtools/samtools \
 && cd samtools \
 && autoreconf -fi \
 && ./configure \
 && make \
 && make install

RUN Rscript -e 'install.packages(c("BiocManager", "remotes", "seqinr", "readr", "tidyr", "dplyr", "snowfall"), repos="https://cloud.r-project.org")' \
 && Rscript -e 'BiocManager::install("Rsamtools", ask=FALSE, update=FALSE)'

RUN git clone https://github.com/KevinLYW366/mlstverse \
 && git clone https://github.com/KevinLYW366/mlstverse.Mycobacterium.db \
 && git clone https://github.com/KevinLYW366/mlstverse.pubmlst.db

RUN Rscript -e 'remotes::install_local("/opt/src/mlstverse", upgrade="never", dependencies=FALSE)' \
 && Rscript -e 'remotes::install_local("/opt/src/mlstverse.Mycobacterium.db", upgrade="never", dependencies=FALSE)' \
 && Rscript -e 'remotes::install_local("/opt/src/mlstverse.pubmlst.db", upgrade="never", dependencies=FALSE)'

RUN rm -rf /opt/src/minimap2 /opt/src/bwa /opt/src/htslib /opt/src/samtools \
           /opt/src/mlstverse /opt/src/mlstverse.Mycobacterium.db /opt/src/mlstverse.pubmlst.db \
 && groupadd -o -g 101 gen-info \
 && useradd -u 10466 -g 101 nanopore-user

WORKDIR /root
ENTRYPOINT ["bash"]
