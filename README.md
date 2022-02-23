recontig
========

[![linux ldc build](https://github.com/blachlylab/recontig/actions/workflows/dbuild-linux.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/dbuild-linux.yml)
[![macos ldc build](https://github.com/blachlylab/recontig/actions/workflows/dbuild-macos.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/dbuild-macos.yml)
[![linux pyd build](https://github.com/blachlylab/recontig/actions/workflows/pybuild-linux.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/pybuild-linux.yml)
[![macos pyd build](https://github.com/blachlylab/recontig/actions/workflows/pybuild-macos.yml/badge.svg)](https://github.com/blachlylab/recontig/actions/workflows/pybuild-macos.yml)
[![codecov](https://codecov.io/gh/blachlylab/recontig/branch/master/graph/badge.svg?token=K8R2FF15EO)](https://codecov.io/gh/blachlylab/recontig)
Due to the divergence in reference naming standards used across popular reference files and genomic databases used in standard bioinformatics analysis, there's was a need for conversion of one database's convention to the other i.e UCSC to NCBI to gencode to ensembl. Use of databases that are siloed in their naming conventions such as dbsnp (NCBI) and gnomad (ensembl) can create downstream issues for workflows and analysis due to contig naming incompatibility. 

*recontig* fufills this need by providing fast conversion of NCBI, UCSC, Ensembl, and Gencode to the other. `recontig` is written in *D* that relies on [dhtslib](https://github.com/blachlylab/dhtslib) - a collection of wrappers and bindings of [htslib](https://github.com/samtools/htslib) for the D language. `recontig` utilizes the power and speed of dhtslib/htslib to convert contig names for `.vcf/.bcf`, `.bam/.sam`, `.gff2/.gff3`, and `.bed` file types. Headers and records are modified properly to specification with support for files compressed with bgzip, gzip, or remote files via http(s) or amazon s3.

### Python - PyD recontig

*Python* support is built through PyD. The speed of *D* is thus kept while enabling the utilites and libraries of *python*. This may be of interest for building `recontig` directly into other projects for quick implementation during tool, workflow development, or analysis. Information about the `recontig` python package can be found [here](INSTALL.md#building-the-recontig-python-package).

### Quick Start
Install `recontig`'s [dependencies](INSTALL.md#dependencies), then download a [binary]() for your system (we support macOS and linux).

Alternatively you can install via bioconda:
```
conda install -c bioconda recontig
```

Or if you would like to build from source, you can find instructions [here](INSTALL.md#building-recontig-from-source).


## Running recontig
`recontig` can read VCF/BCF, SAM/BAM, GFF, BED, and custom formats that are delimited record-based text formats. All outputs are written to stdout (unless using the `-o` flag) as the uncompressed, text based format of the original file type. The only exception to this is using `-f sam` will output SAM while using using `-f bam` will output BAM. All file inputs can be remote (accessed via web link), gzipped, or bgzipped. `recontig` will handle downloading and decompression.  
### General usage:
Usage of a dpryan79 mapping file (automatically downloaded on-the-fly)
```
# output is printed to stdout 
./recontig convert -b GRCh37 -c UCSC2ensembl -f vcf in.vcf.gz > out.vcf
```
Usage of a specific mapping file
```
./recontig convert -m mapping.txt.gz -f bed in.bed > out.bed
```
Usage of a generic file format
```
# defaults if not supplied --delimiter '\t' --comment '#'
./recontig convert -m mapping.txt.gz --col 1 --delimiter ',' --comment '#' in.txt > out.txt
```

Usage of a SAM/BAM file outputing BAM
```
./recontig convert -m mapping.txt -f bam in.bam > out.bam
```

Usage of a SAM/BAM file outputing SAM
```
./recontig convert -m mapping.txt -f sam in.bam > out.sam
```

Web-based access (try me)
```
./recontig convert -m https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_ensembl2UCSC.txt -f vcf https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.Y.vcf.bgz | less -S
```

**Note**: Due to the nature of contig renaming, all resultant files from `recontig` will need to be resorted using the appropriate tool (i.e `samtools sort`, `bcftools sort`, `bedtools sort`).
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
./recontig conversion-help -b selected-build
```

### CLI
```
recontig: convert contig names for different bioinformatics file types.

Subcommands:
build-help          check availiable builds from dpryan79's github
conversion-help     check availiable conversions for a specified build from dpryan79's github
convert             convert a file from one naming convention to another
make-mapping        make a contig conversion file from two fasta files
```
```
recontig conversion-help: check availiable conversions for a specified build from dpryan79's github

usage: recontig conversion-help -b build 

-b   --build Genome build i.e GRCh37 for using dpryan79's files
-q   --quiet silence warnings
-v --verbose print extra information
     --debug print extra debug information
-h    --help This help information.
```
```
recontig convert: remap contig names for different bioinformatics file types.

usage: recontig convert [-e ejected.txt] [-o output] [-m mapping.txt | -b build -c conversion] [-f filetype | --col 1 --delimiter ','] <in.file>

Input can be any of the following formats: vcf, bcf, bam, sam, bed, gff
Input can also be a delimited record based file 
Input can be compressed with gzip or bgzf and can be accessed remotely via https or s3 (see htslib for details).

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
```
recontig make-mapping: make a contig conversion file from two fasta files
Makes a mapping file for two fasta files that have been faidx'd (samtools faidx)
Fastas can be compressed with bgzf and can be accessed remotely via https or s3 (see htslib for details).

usage: recontig make-mapping [-o output] <from.fa> <to.fa>

-o  --output name of file out (default is - for stdout)
-q   --quiet silence warnings
-v --verbose print extra information
     --debug print extra debug information
-h    --help This help information.
```

## Common Problems and solutions
### Python versions
On macOS it seems `python3` works best. 
