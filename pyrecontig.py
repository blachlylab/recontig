import pandas as pd
import recontig
import argparse
import sys
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
 
def lengthCheck(pd1, pd2, gft):
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
    print("Checking length of mappings for errors", file = sys.stderr)
    # Sorting the information from the loaded gft
    chroms = gft.iloc[:,0]
    starts = gft.iloc[:,3]
    stops = gft.iloc[:,4]

    pair = []
    buildDict = {}
    # Add  the keys to the dictionary
    buildDict = {k:[] for k in chroms}
    # Add the largest range for a start and stop for a componenet as defined by ucsc
    for i in range(0,len(chroms)):
        pair = buildDict[chroms[i]]
        # add only the minimum and maxmimum to the pairing
        if len(pair) == 0:
            pair.append(starts[i])
            pair.append(stops[i])
        else:
            if starts[i] < pair[0]:
                pair[0] = starts[i]
            if stops[i] > pair[1]:
                pair[1] = stops[i]
    
    # dattern for chromosomes to set the index to 0 instead of 1000 as the start
    pattern = ["chr" + str(x) for x in range(1,23)]
    pattern.append("chrX")
    pattern.append("chrY")
    pattern.append("chrM")

    # Correct for the front end length overhang.
    for key in pattern:
        try:
            pair = buildDict[key]
            pair[0] = 0
            buildDict[key] = pair
        except:
            print("Warning: " + key  + " not found in hg38 ucsc reference.", 
                    file = sys.stderr)
    
    # Check for chromosome discrpecnies
    for key in buildDict.keys():
        pair = buildDict[key]
        chromInDictPd = pd2[pd2["#CHROM"].str.contains(key)]
        # Any value that is not found is simply not found in the input vcf from the user 
        # given for conversion.
        if not len(chromInDictPd) == 0:
            pos = list(chromInDictPd["POS"])
            sizeOut = [x for x in pos if int(x) > int(pair[1])]
            # Output warning for those not fit between interval.
            if not len(sizeOut) == 0:
                print("Warning: Unable to varify: " + str(len(sizeOut)) + " out of " + str(len(pos)) + " on component " + key, 
                        file = sys.stderr)
    print("Done validating chromosome lengths", file = sys.stderr)
    
    outLST.append(pd1)
    outLST.append(pd2)

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
