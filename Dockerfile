FROM ubuntu:20.04

LABEL description="Docker image for recontig with PyD"

ARG DEBIAN_FRONTEND=noninteractive
# Change htslib version here
ENV htslib_ver 1.13  
ENV PACKAGES autoconf automake make gcc perl pip git \
    curl zlib1g-dev libbz2-dev liblzma-dev libssl-dev \
    libcurl4-gnutls-dev libxml2

# Install htslib Dependencies and library
# Install packages
RUN apt-get update && apt-get install \
    && apt-get install -y ${PACKAGES} && apt-get clean

ADD https://github.com/samtools/htslib/releases/download/${htslib_ver}/htslib-${htslib_ver}.tar.bz2 /usr/local/
RUN cd /usr/local/ && tar -xvf /usr/local/htslib-${htslib_ver}.tar.bz2 && \
    cd htslib-${htslib_ver} && ./configure && make && make install

# Install the ldc
RUN curl https://dlang.org/install.sh | bash -s ldc

# Install CPython and pyd
RUN pip install cython pyd

# Install recontig
RUN cd usr/local/ && git clone --recurse-submodules https://github.com/blachlylab/recontig.git
RUN . ~/dlang/ldc-*/activate && cd /usr/local/recontig \
    && LD_LIBRARY_PATH=/usr/local/htslib-${htslib_ver} LIBRARY_PATH=/usr/local/htslib-${htslib_ver} dub build

# Define default location
ENTRYPOINT ["/usr/local/recontig/recontig"]