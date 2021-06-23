# Install and Setup
## Dependencies
### htslib
`recontig` and by extension [dhtslib](https://github.com/blachlylab/dhtslib) rely on [htslib](https://github.com/samtools/htslib). htslib relies on a handful of compression and web-access libraries: zlib, bzip2, lzma, curl, and ssl. Technically htslib can be built without some of these libraries, though to get all the benefits of `recontig` we recommend installing all of them.

To intall htslib dependencies:
```
Debian / Ubuntu
---------------

sudo apt-get update  # Ensure the package list is up to date
sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev

Note: libcurl4-openssl-dev can be used as an alternative to libcurl4-gnutls-dev.

RedHat / CentOS / Amazon Linux
---------------

sudo yum install autoconf automake make gcc perl-Data-Dumper zlib-devel bzip2 bzip2-devel xz-devel curl-devel openssl-devel

Alpine Linux
------------

sudo apk update  # Ensure the package list is up to date
sudo apk add autoconf automake make gcc musl-dev perl bash zlib-dev bzip2-dev xz-dev curl-dev libressl-dev

OpenSUSE
--------

sudo zypper install autoconf automake make gcc perl zlib-devel libbz2-devel xz-devel libcurl-devel libopenssl-devel

MacOS
-----

brew install xz autoconf automake
```

You will then need to download [htslib](http://www.htslib.org/download/).
As of now we support versions >=1.10. To install htslib:
```
curl https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 -o htslib-1.12.tar.bz2
tar -xjf htslib-1.12.tar.bz2
cd htslib-1.12
./configure
make 
sudo make install
```

## Building recontig from source
Clone the repository via 
```
git clone --recurse-submodules https://github.com/blachlylab/recontig.git
```
If you already have htslib 1.12 installed on your system and on your path then activate your *D* environment
```
dub build
```
If your htslib install is in a non-standard location
```
LD_LIBRARY_PATH=path/to/htslib/ LIBRARY_PATH=path/to/htslib/ dub build
```

## Building the recontig python package
### Cython and PyD
Install `CPython >=3.7` [here](https://www.python.org/downloads/) is it's not already on your system. [PyD](https://github.com/ariovistus/pyd) does not have shared library support for GDC, so it is recommended to install using `DMD` or `LDC fe2.065+`
Best way to get both:
```
pip install cython pyd

#macOS
pip3 install cython pyd
```

To install *D* run the following command to install into your home folder. 
```
curl https://dlang.org/install.sh | bash -s ldc
```
To activate the *D* environment: `source ~/dlang/ldc-*/activate`

You could use `dmd` instead of `ldc`.
```
curl https://dlang.org/install.sh | bash -s dmd
```
We prefer `ldc` for its better performance and it is the compiler we actively test `recontig` with.
Your results may vary.

### setup.py build
If not, then we have provided an install script that will install htslib on your system or check if it was provided if the user has proper system permissions to do so. We recommend running this with the `ldc`. Extended information about the python package can be found [here]().
```
python setup.py build -c ldc
python setup.py install
#macos
python3 setup.py build -c ldc
python3 setup.py install
```
To test your install:
```
python -c "import recontig"
#macos
python3 -c "import recontig"
```