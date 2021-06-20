import os
import urllib.request
import pandas as pd
import recontig
import argparse
import sys
import gzip


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

def _ungzip(vcf):
    """
        Takes in a gziped file and opens it assuming it ends in '.gz'.
        The output is then an uncompressed file reading for reading by
        other functions.
        Input: path to vcf file as string
        Output: a decompressed file
    """
    # Open the file with the gzip package.
    if vcf.endswith(".gz"):
        # For a vcf that is gzipped
        vcfDecomp = gzip.open(vcf, 'rt')
        try:
            vcfDecomp.read(1)
            print("Uncompressed vcf file: " + vcf, file=sys.stderr)
            return vcfDecomp
        except gzip.BadGzipFile:
            print("GZIPFile is a bad zip file: " + vcf, file=sys.stderr)
    else:
        # for a vcf that is not compressed
        vcf = open(vcf, 'rt')
        return vcf

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


def _lengthCheckVcf(pd1, pd2, gft):
    """
       Uses a gtf file from each respective build and version to determine 
       the lengths. These are then used as a check for positions that fall out
       of bounds during the switch. In particular, this is used to check the
       mitochondion since the build will vary depending on the version and
       may need to be adjusted.
       Input: Two pandas frames with pd1 containing the header vcf and pd2 
            containing the variants.
       Input: url linking to the gtf file to download.
       Output: Backfilled pandas frames in a list with position one containing
            pd1 and position two containing pd2.
    """
    outLST = []
    print("Checking length of converted contigs for errors", file = sys.stderr)
    # Sorting the information from the loaded gft
    chroms = gft.iloc[:,0]
    starts = gft.iloc[:,3]
    stops = gft.iloc[:,4]

    pair = []
    buildDict = {}

    # Add  the keys to the dictionary
    buildDict = {k:[] for k in chroms}

    # Add the largest range for a start and stop for a componenet as defined by 
    # given gtf. Nothing thus will fall outside of this range assuming the 
    # chromosomes have been sucessfully converted over.
    for i in range(0,len(chroms)):
        pair = buildDict[chroms[i]]

        # Add only the minimum and maxmimum to the pairing
        if len(pair) == 0:
            pair.append(starts[i])
            pair.append(stops[i])
        else:
            if starts[i] < pair[0]:
                pair[0] = starts[i]
            if stops[i] > pair[1]:
                pair[1] = stops[i]
    print("...sucessfully built max-min dictionary", file = sys.stderr)

    # Sort on chromosome and position columns then reset the index.
    pd2['POS'] = pd2['POS'].astype('float64')
    convertedMin = pd2.loc[pd2.groupby('#CHROM')['POS'].idxmin()]
    convertedMin.rename(columns={'POS':'MIN'}, inplace=True)
    convertedMax = pd2.loc[pd2.groupby('#CHROM')['POS'].idxmax()]
    convertedMax.rename(columns={'POS':'MAX'}, inplace=True)   
    print("...sucessfully built max-min for each contig", file = sys.stderr)

    # Check against the unchanged lengths of the chromosomes in the header
    contigPd = pd1[pd1['head'].str.contains('##contig=')]
    contigPd = contigPd['head'].str.split(',', expand=True)
    chromPd = contigPd.iloc[:, 0].str.split("##contig=<ID=", expand=True)
    lenPd = contigPd.iloc[:, 1].str.split('=', expand=True)
    contigPd["#CHROM"] = chromPd.iloc[:, 1]
    contigPd["#LEN"] = lenPd.iloc[:, 1] 
    contigPd = contigPd[["#CHROM", "#LEN"]]
    print("...collected and sorted contig information from file", file = sys.stderr)

    # Sort on the chromosomes to ensure same order for comparison
    contigPd = contigPd.sort_values('#CHROM', ascending=True)
    convertedMax = convertedMax.sort_values('#CHROM', ascending=True)

    # The length list for all contigs, this is unchanged from conversion
    # in the D package, thus is used as a check if length is the same 
    # for the contigs that have been converted.
    contigPd = contigPd.iloc[:, 1] 
    contigsLST = convertedMax['#CHROM'].tolist()
    lengthLST = contigPd.tolist()
    convertedMaxLST = convertedMax["MAX"].tolist()
    print("...generated length lists", file = sys.stderr)

    if not len(lengthLST) == len(convertedMaxLST):
        print("WARNING: list lengths are not same size in _lengthCheckVcf()", file = sys.stderr)
        return False
    
    for i in range(0,len(lengthLST)):
        if int(lengthLST[i]) < int(convertedMaxLST[i]):
            print("WARNING: contig " + contigsLST[i] + " does not match position that is expected in the coverted vcf")
            print("position before conversion: " + str(lengthLST[i]) + ". After conversion: " + str(convertedMaxLST))
            print("Coversion contig falls out of original length by " + str(abs(int(lengthLST[i])-int(convertedMaxLST))))
    print("...length check done!", file = sys.stderr)
    
    return True


def main():
    # GRCh38 GTF Files.
    ucscGft = "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.refGene.gtf.gz" 
    gencodeGft = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz"
    ensemblGft = "http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
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
    gft = {} # Empty initialization
    if(args.mapping != None):
        mapping = _getmapping(args.mapping)
    elif(args.build and args.conversion):
        mapping = _getdpyryan(args.build, args.conversion)
        # Download correct gtf files from the right most argument.
        whichGtf = args.conversion.split('2')[1]
        if "UCSC" in whichGtf:
            gft = _gftRead(ucscGft, 0)
        elif "gencode" in whichGtf:
            gft = _gftRead(gencodeGft, 0)
        elif "ensembl" in whichGtf:
            gft = _gftRead(ensemblGft, 5)
        elif "refSeq" in whichGtf:
            gft = _gftRead(refseqGft, 0)
        print("Downloaded GTF needed for length checks", file = sys.stderr)
    else:
        raise Exception("Please provide either mapping file or build and conversion")

    # For given argument, run the conversion.
    if args.fileType == "vcf":
        # Convert the vcf over to the desired naming convention.
        recontig.recontigVcf(args.file,"ejected.vcf", mapping, name, "")
        #convertedVcf = _ungzip(args.output)
        convertedVcf = open(args.output, 'r')
        vcfFrame = _vcfReadToPandas(convertedVcf)
        # Check lengths of vcf after conversion.
        if _lengthCheckVcf(vcfFrame[0], vcfFrame[1], gft) == True:
            print("Length check passed.", file = sys.stderr)
        else:
            print("WARNING: Length check failed.", file = sys.stderr)

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
