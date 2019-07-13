FROM rnakato/ubuntu:18.04
LABEL original from duplexa/4dn-hic, modified by Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>

# 1. general updates & installing necessary Linux components
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    bzip2 \
    default-jdk \
    gawk \
    gcc \
    git \
    less \
    liblz4-tool \
    libncurses-dev \
    libxml2-dev \
    make \
    python3.6-dev \
    python3-setuptools \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://bootstrap.pypa.io/get-pip.py --no-check-certificate \
    && python3.6 get-pip.py

# installing R & dependencies for pairsqc
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    r-base \
    r-base-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R CMD javareconf \
    && R -e 'install.packages("devtools", repos="https://cran.ism.ac.jp/")' \
    && R -e 'install.packages( "Nozzle.R1", type="source", repos="https://cran.ism.ac.jp/" )' \
    && R -e 'library(devtools); install_url("https://github.com/SooLee/plotosaurus/archive/0.9.2.zip")' \
    && R -e 'install.packages("stringr", repos="https://cran.ism.ac.jp/" )'

# installing conda
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && bash Miniconda2-latest-Linux-x86_64.sh -p /miniconda2 -b
ENV PATH=/miniconda2/bin:$PATH
RUN conda update -y conda \
    && rm Miniconda2-latest-Linux-x86_64.sh

# installing gawk for juicer
RUN echo 'alias awk=gawk' >> ~/.bashrc

# download tools
WORKDIR /usr/local/bin
COPY downloads.sh .
RUN . downloads.sh

# set path
ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH
ENV PATH=/usr/local/bin/pairix/bin/:/usr/local/bin/pairix/util/:$PATH
ENV PATH=/usr/local/bin/pairix/util/bam2pairs/:$PATH
ENV PATH=/usr/local/bin/pairsqc/:$PATH
ENV PATH=/usr/local/bin/juicer/CPU/:/usr/local/bin/juicer/CPU/common:$PATH
ENV PATH=/usr/local/bin/hic2cool/:$PATH
ENV PATH=/usr/local/bin/mcool2hic/:$PATH
ENV PATH=/usr/local/bin/FastQC/:$PATH

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# wrapper
COPY scripts/ .
RUN chmod +x run*.sh

# default command
CMD ["run-list.sh"]