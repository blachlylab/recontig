# Recontig
Due to the divergence in reference naming standards in commonly needed file used in standard bioinformatics analysis, there's was a need for conversion of one database's convention to the other i.e UCSC to NCBI to gencode to ensembl. Use of databases that are siloed in their naming conventions such as dbsnp (NCBI) and gnomad (ensembl) can create downstream issues for workflows and analysis due to contig naming incompatibility. 

*Recontig* fufills this need by providing fast conversion of NCBI, UCSC, Ensembl, and Gencode to the other. Recontig is written in *D* that relies on [dhtslib](https://github.com/blachlylab/dhtslib) - a collection of wrappers and bindings of htslib for the D language. Recontig utilizes the power and speed of dhtslib/htslib to convert contig names for `.vcf/.bcf`, `.bam/.sam`, `.gff2/.gff3`, and `.bed` file types. Headers and records are modified properly to specification with support for files compressed with bgzip, gzip, or remote files via http(s) or amazon s3.

### Python - PyD Recontig

*Python* support is built through PyD. The speed of *D* is thus kept while enabling the utilites and libraries of *python*. This may be of interest for building recontig directly into other projects for quick implementation during tool, workflow development, or analysis. 

## Install and Setup
### Dependencies
Install `CPython >=3.7` [here](https://www.python.org/downloads/) is it's not already on your system.

[PyD](https://github.com/ariovistus/pyd) does not have shared library support for GDC, so it is recommended to install using `DMD` or `LDC fe2.065+`

You will also need you install `wget` if it is not already. `pip3 install wget` should do this for you.

To install *D* run the following command where every you with the dlang directory to reside. 
```
mkdir -p ~/dlang && wget https://dlang.org/install.sh -O ~/dlang/install.sh
```
To activate if it is not already then type `~/dlang/install.sh install <compiler> -a`

[dhtslib](https://github.com/blachlylab/dhtslib) must be installed as well as [htslilb](http://www.htslib.org/download/).

Finally, PyD can be installed via pip `pip install pyd`.
### Installing on Linux - Ubuntu
Clone the repository via 
```
git clone --recurse-submodules https://github.com/blachlylab/recontig.git
```
If you already have htslib 1.12 installed on your system and on your path then you can run the following for good measure to build
```
LD_LIBRARY_PATH=$HOME/tools/htslib/ LIBRARY_PATH=$HOME/tools/htslib/ dub build
```
If not, then we have provided an install script that will install htslib on your system or check if it was provided if the user has proper system permissions to do so. We recommend running this with Dlang-ldc. 
```
python setup.py build -c ldc
```
## Running recontig

## Common Problems and solutions
### Python versions
Sometimes python 2 will be the default version to run on a system. This may be the case if you are recieving an error like:
```
Traceback (most recent call last):
  File "setup.py", line 110, in <module>
    dhtslibSources = glob.glob(os.path.join("dhtslib","source","**","*.d"), recursive=True)
TypeError: glob() got an unexpected keyword argument 'recursive'
```
when running the `setup.py` as Recontig requires python 3.7 and greater. To solve this consider installing python 3 and setting to an alias like `alias python=python3`
