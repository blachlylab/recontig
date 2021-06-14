import pandas as pd
import recontig
import argparse
import sys
import os
import wget

def _getUserArgs():
    """ Collect user arguments and return an
    argparse object.
    Input: User input through argparse
    Output: Argparse object
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--build","-b",type=str,
            help="Build of file, i.e grch38 or grch37")
    parser.add_argument("--conversion","-c",type=str,
            help="which organization naming convention to convert to")
    parser.add_argument("--fileType","-t",type=str,
            help="The type of file to convert - vcf, gff, bed, bam")
    parser.add_argument("--file", "-f",type=str,
            help="Input file")
    args = parser.parse_args()   

    return args

def _getdpyryan(build, conversion):
    """ Download dpryan's mapping files for the
    quality cross-checking.
    Input: build - the build the file current is written in
    Input: conversion - what build to convert to
    Output: File for dpryan conersion
    """
    return recontig.getdpyryan(build, conversion)

def _getmapping(dpyryanMap):
    """ Read mapping file into pandas frame from the 
    Dlang function 'getmapping'. 
    Input: string to dpyryanMap downloaded file
    Output: mapping file in a pandas' frame
    """
    mapping=pd.read_csv(recontig.getmapping(dpyryanMap), sep="\t")
    return mapping

def _gftRead(url, step):
    """
        Reads in a gtf file from a specific db given the url.
        Some gft have a certain number of header lines that 
        are skipped however.
        Input: url where gtf is fetched from
        Input: number of lines to skip while reading in the frame
        Output: gtf in a pandas frame
    """
    vcf = wget.download(url, out = "/tmp/conversion.gtf.gz")
    gtf = pd.read_csv(vcf,
            compression = "gzip",
            engine = "python",
            delimiter = '\t',
            skiprows = step,
            header = None)

    return gtf

def _createFileOutName(fileName, ext):
    """ Creates the default output file name from the given file
    Input: filename - string representing path and filename
    Input: ext - extension to split
    Output: Outputfile on same path as input
    """
    name = os.path.basename(fileName)
    dirname = os.path.dirname(fileName)

    strippedName = name.split('.')
    lastE = len(strippedName) # Last element
    # Include all variables but .ext and .ext.gz
    if strippedName[lastE-1] == ext:
        nameLST = strippedName[0:len(strippedName)-1]
        nameLST.append("converted." + ext)
        name =  "".join(nameLST)
    elif strippedName[lastE-2] == ext:
        nameLST = strippedName[0:len(strippedName)-2]
        nameLST.append("converted." + ext + ".gz")
        name =  ".".join(nameLST)        
    # Connect directory and new output name.
    out = os.path.join(dirname, name)

    return out
    

def main():
    # GRCh38 GTF Files.
    ucscGft = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.refGene.gtf.gz" 
    gencodeGft = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz"
    ensemblGft = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz"
    refseq = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
    
    # Get arguments from user.
    args = _getUserArgs()
    name = _createFileOutName(args.file, args.fileType)
    # Get dpyryans files for cross-checks.
    dpyryan = _getdpyryan(args.build, args.conversion)
    # Get mapping
    mapping = _getmapping(dpyryan)   
    # Download correct gtf files
    gft = {} # Empty initialization
    whichGtf = args.conversion.split('2')
    if "UCSC" in whichGtf:
        gft = _gftRead(ucscGft, 0)
    elif "GENCODE" in whichGtf:
        gft = _gftRead(gencodeGft, 0)
    elif "ENSEMBL" in whichGtf:
        gft = _gftRead(ensemblGft, 0)
    elif "REFSEQ" in whichGtf:
        gft = _gftRead(refseqGft, 0)

    # For given argument, run the conversion.
    if args.fileType == "vcf":
        recontig.recontigVcf(args.filename,name,{"1":"chr1"},"helpme")
    elif args.fileType == "bed":
        recontig.recontigBed(args.filename,name,{"1":"chr1"},"helpme")
    elif args.fileType == "bam":
        recontig.recontigBam(args.filename,name,{"1":"chr1"},"helpme")
    elif args.fileType == "gff":
        recontig.recontigGff(args.filename,name,{"1":"chr1"},"helpme")


if __name__ == "__main__":
    main()
