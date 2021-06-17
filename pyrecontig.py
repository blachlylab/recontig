import pandas as pd
import recontig
import argparse
import csv
import urllib.request


def _getdpyryan(build, conversion):
    """ Download dpryan's mapping files for the
    quality cross-checking.
    Input: build - the build the file current is written in
    Input: conversion - what build to convert to
    Output: mapping dictionary
    """
    return recontig.getDpryan79ContigMapping(build, conversion)

def _getmapping(mappingFile):
    """ Read local mapping file into dictionary. 
    Input: string to local mapping file
    Output: mapping dictionary
    """
    mapping = recontig.getContigMapping(mappingFile)
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
    urllib.request.urlretrieve(url, "/tmp/conversion.gtf.gz")
    gtf = pd.read_csv("/tmp/conversion.gtf.gz",
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
 
def _vcfReadToPandas(vcf):
    """ Takes a vcf file and reads it into a pandas data frame.
       Input: vcf file.
       Output: A list containing two pandas frames. The first containing the
            vcf header and the second the data.
    """
    dataDict = {}
    finHeadLST = []; finDataLST = []; outLST = []
    # The vcf is separated into different parts, with the header in one frame 
    # and the variants in another.
    for line in vcf:
        # Add the lines for the vcf header to a separate frame by building
        # them into dictionaries + lists. After they will be converted below
        # to the pandas data frame. The same occurs for the data portion of 
        # the vcf file. 
        if line.startswith("##"):
            headDict = {}
            # Head is used as the name of the header frame column.
            headDict["head"] = line.strip()
            finHeadLST.append(headDict)
        if not line.startswith("##"):
            header = []
            if line.startswith("#") and not line.startswith("##"):
                # Create the header for just the variant frame.
                header = line.strip().split('\t')
                dataDict = {key:[] for key in header}
            else:
                # All other values are used for the variant frame section.
                line = line.strip().split('\t')
                keys = list(dataDict.keys())
                # In python 3.7 upwards dictionaries maintain order. Here we
                # fill this in order of lines with retenion to order only noted
                # to keep an ordered vcf in the output. In other versions of python
                # an unordered vcf will be output.
                for cnt in range(0,len(keys)):
                    dataDict[keys[cnt]] = line[cnt]
                    if cnt == len(keys)-1:
                        # Break the duplicte referencing structure by copying 
                        # the dictionary.
                        finDataLST.append(dataDict.copy())

    # Set the head frame.
    headFrame = pd.DataFrame(finHeadLST, columns=["head"])
    # Set the data frame with specified columns to preserve the order.
    dataFrame = pd.DataFrame(finDataLST)
    
    outLST.append(headFrame)
    outLST.append(dataFrame)
    
    return outLST


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
    parser.add_argument("--fileType","-t",type=str, required=True,
            help="The type of file to convert - vcf, gff, bed, bam")
    parser.add_argument("--file", "-f",type=str, required=True,
            help="Input file")
    parser.add_argument("--mapping", "-m",type=str,
        help="Name of mapping file")
    parser.add_argument("--output", "-o",type=str,
            help="ouput file", default="-")
    args = parser.parse_args()

    return args


def writeOutVcf(vcfFrame, out):
    """ Writes the vcfFrame into a vcf file.
        Input: List containing two frames with the
            vcf header in one frame and the variant 
            data in the other.
        Input: Output file to be written to.
    """
    # Write out the headers
    vcfFrame[0].to_csv(out,
            index = False,
            header = False,
            quoting = csv.QUOTE_NONE,
            escapechar = '\t',
            sep = '\t',
            line_terminator = '\n')
    # Write out the variants
    vcfFrame[1].to_csv(out,
            sep = '\t',
            index = False,
            header = True)

    out.close()


def main():
    # GRCh38 GTF Files.
    ucscGft = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.refGene.gtf.gz" 
    gencodeGft = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz"
    ensemblGft = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz"
    refseqGft = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz"
    
    # Get arguments from user.
    args = _getUserArgs()

    # Get the outfie name if it was not supplied.
    name = args.output
    if args.output is None:
        name = _createFileOutName(args.file, args.fileType)

    # Get dpyryans files for cross-checks.
    # Get mapping
    mapping = {}
    if(args.mapping != None):
        mapping = _getmapping(args.mapping)
    elif(args.build and args.conversion):
        mapping = _getdpyryan(args.build, args.conversion)
        # Download correct gtf files
        gft = {} # Empty initialization
        whichGtf = args.conversion.split('2')[0]
        if "UCSC" in whichGtf:
             gft = _gftRead(ucscGft, 0)
        elif "gencode" in whichGtf:
             gft = _gftRead(gencodeGft, 0)
        elif "ensembl" in whichGtf:
             gft = _gftRead(ensemblGft, 0)
        elif "RefSeq" in whichGtf:
             gft = _gftRead(refseqGft, 0)
    else:
        raise Exception("Please provide either mapping file or build and conversion")

    # For given argument, run the conversion.
    if args.fileType == "vcf":
        recontig.recontigVcf(args.file,"ejected.vcf", mapping, name, "")
        convertedVcf = open(args.output, 'r')
        vcfFrame = _vcfReadToPandas(convertedVcf)
        writeOutVcf(vcfFrame, open(name, "w"))
    elif args.fileType == "bed":
        recontig.recontigBed(args.file,"ejected.bed",mapping, args.output, "")
    elif args.fileType == "bam":
        recontig.recontigBam(args.file,"ejected.sam",mapping, args.output, "")
    elif args.fileType == "sam":
        recontig.recontigSam(args.file,"ejected.sam",mapping, args.output, "")
    elif args.fileType == "gff":
        recontig.recontigGff(args.file,"ejected.gff",mapping, args.output, "")

if __name__ == "__main__":
    main()
