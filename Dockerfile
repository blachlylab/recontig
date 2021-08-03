FROM ubuntu:20.04

# Install htslib Dependencies and library
RUN apt-get update
RUN apt-get install -y autoconf automake make gcc \ 
    perl zlib1g-dev libbz2-dev liblzma-dev git \
    libcurl4-gnutls-dev libssl-dev wget curl pip libxml2
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
RUN tar -xvf htslib-1.13.tar.bz2 && cd htslib-1.13 && ./configure && make && make install

# Add To Path
RUN LD_LIBRARY_PATH=$PWD/htslib-1.13 LIBRARY_PATH=$PWD/htslib-1.13

# Install CPython
RUN pip install cython pyd
# Install the ldc
RUN curl https://dlang.org/install.sh | bash -s ldc
# Clone recontig
RUN git clone --recurse-submodules https://github.com/blachlylab/recontig.git
# Compile recontig

#RUN . ~/dlang/ldc-1.27.0/activate && cd recontig \
#    && python3 setup.py build -c ldc && setup.py install