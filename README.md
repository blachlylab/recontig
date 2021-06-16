recontig
========

[![linux ldc build](https://github.com/blachlylab/recontig/actions/workflows/dbuild-linux.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/dbuild-linux.yml)
[![macos ldc build](https://github.com/blachlylab/recontig/actions/workflows/dbuild-macos.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/dbuild-macos.yml)
[![linux pyd build](https://github.com/blachlylab/recontig/actions/workflows/pybuild-linux.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/pybuild-linux.yml)
[![macos pyd build](https://github.com/blachlylab/recontig/actions/workflows/pybuild-macos.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/pybuild-macos.yml)

Due to the divergence in reference naming standards used across popular reference files and genomic databases used in standard bioinformatics analysis, there's was a need for conversion of one database's convention to the other i.e UCSC to NCBI to gencode to ensembl. Use of databases that are siloed in their naming conventions such as dbsnp (NCBI) and gnomad (ensembl) can create downstream issues for workflows and analysis due to contig naming incompatibility. 

*recontig* fufills this need by providing fast conversion of NCBI, UCSC, Ensembl, and Gencode to the other. `recontig` is written in *D* that relies on [dhtslib](https://github.com/blachlylab/dhtslib) - a collection of wrappers and bindings of [htslib](https://github.com/samtools/htslib) for the D language. `recontig` utilizes the power and speed of dhtslib/htslib to convert contig names for `.vcf/.bcf`, `.bam/.sam`, `.gff2/.gff3`, and `.bed` file types. Headers and records are modified properly to specification with support for files compressed with bgzip, gzip, or remote files via http(s) or amazon s3.

### Python - PyD recontig

*Python* support is built through PyD. The speed of *D* is thus kept while enabling the utilites and libraries of *python*. This may be of interest for building `recontig` directly into other projects for quick implementation during tool, workflow development, or analysis. 

## Install and Setup
### Dependencies
#### htslib
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

#### Cython and PyD
Install `CPython >=3.7` [here](https://www.python.org/downloads/) is it's not already on your system. [PyD](https://github.com/ariovistus/pyd) does not have shared library support for GDC, so it is recommended to install using `DMD` or `LDC fe2.065+`
Best way to get both:
```
pip install cython pyd

#macOS
pip3 install cython pyd
```

To install *D* run the following command where every you with the dlang directory to reside. 
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

### Building recontig
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
### Building the recontig python package
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
## Running recontig
`recontig` can read VCF/BCF, SAM/BAM, GFF, BED, and custom formats that are delimited record-based text formats. All outputs are written to stdout (unless using the `-o` flag) as the uncompressed, text based format of the original file type. The only exception to this is using `-f sam` will output SAM while using using `-f bam` will output BAM. All file inputs can be remote (accessed via web link), gzipped, or bgzipped. `recontig` will handle downloading and decompression.  
### General usage:
Usage of a dpryan79 mapping file (automatically downloaded on-the-fly)
```
# output is printed to stdout 
./recontig -b GRCh37 -c UCSC2ensembl -f vcf in.vcf.gz > out.vcf
```
Usage of a specific mapping file
```
./recontig -m mapping.txt.gz -f bed in.bed > out.bed
```
Usage of a generic file format
```
# defaults if not supplied --delimiter '\t' --comment '#'
./recontig -m mapping.txt.gz --col 1 --delimiter ',' --comment '#' in.txt > out.txt
```

Usage of a SAM/BAM file outputing BAM
```
./recontig -m mapping.txt -f bam in.bam > out.bam
```

Usage of a SAM/BAM file outputing SAM
```
./recontig -m mapping.txt -f sam in.bam > out.sam
```

Web-based access (try me)
```
./recontig -m https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_ensembl2UCSC.txt -f vcf https://storagegoogleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz | less -S
```
### Make a mapping file
`recontig` can create mapping files by comparing two faidx'd fasta files. All contigs are compared for matching md5sums. Contigs with matching sums are reported in the output. This output can then be used with `recontig` to convert supported files.
```
samtools faidx UCSC.fasta #could be bgzipped
samtools faidx ensembl.fasta
./recontig make-mapping UCSC.fasta ensembl.fasta > UCSC2ensembl.txt
```

### Check build and conversion options
`recontig` downloads files from dpryan79's [ChromosomeMappings](https://github.com/dpryan79/ChromosomeMappings) github repository.
To check the availiable builds that are availiable:
```
./recontig build-help
```
To check the availiable conversions for a particular build that are availiable:
```
./recontig -b selected-build conversion-help
```

### CLI
```
recontig: remap contig names for different bioinformatics file types.

usage: recontig [-e ejected.txt] [-o output] [-m mapping.txt | -b build -c conversion] [-f filetype | --col 1 --delimiter ','] <in.file>

Input can be any of the following formats: vcf, bcf, bam, sam, bed, gff
Input can also be a delimited record based file
Input can be compressed with gzip or bgzf and can be accessed remotely via https or s3 (see htslib for details).
use 'recontig build-help' to check availiable builds
use 'recontig -b build' conversion-help to check availiable conversions for a build
use 'recontig make-mapping <from.fasta> <to.fasta>' to make a mapping file from two faidx'd fasta files

-b          --build Genome build i.e GRCh37 for using dpryan79's files
-c     --conversion Conversion string i.e UCSC2ensembl for using dpryan79's files
-e --ejected-output File to write ejected records to (records with unmapped contigs)
-f      --file-type Type of file to convert (vcf, bcf, bam, sam, bed, gff)
-m        --mapping If want to use your own remapping file instead of dpryan79's
-o         --output name of file out (default is - for stdout)
-q          --quiet silence warnings
-v        --verbose print extra information
              --col if converting a generic file you can specify a column
          --comment if converting a generic file you can specify what a comment line starts with (default: '#')
            --debug print extra debug information
        --delimiter if converting a generic file you can specify a delimiter (default: '\t')
-h           --help This help information.
```


## Common Problems and solutions
### Python versions
On macOS it seems `python3` works best. 
