#from recontigPyD import vcf, bed, gff, bam
import recontig
import argparse
import sys

def getUserArgs():
    """ Collect user arguments and return an
    argparse object.
    Input: User input through argparse
    Output: Argparse object
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--build","-b",type=str)
    parser.add_argument("--conversion","-c",type=str)
    parser.add_argument("filename")
    args = parser.parse_args()   

    return args

def getdpyryan(build, conversion):
	return recontig.getdpyryan(build, conversion)

def getmapping(dpyryanMap):
	return recontig.getmapping(dpyryanMap)

def main():
    # Get arguments from user.
    args = getUserArgs()
    # Get dpyryans files for cross-checks.
    dpyryan = getdpyryan(args.build, args.conversion)
    # Get mapping
    mapping = getmapping(dpyryan)
    # For given argument, run the conversion.
    if args.filename == "vcf":
        recontig.recontigVcf(args.filename,"ejected.vcf",{"1":"chr1"},"helpme")
    elif args.filename == "bed":
        recontig.recontigBed(args.filename,"ejected.bed",{"1":"chr1"},"helpme")
    elif args.filename == "bam":
        recontig.recontigBam(args.filename,"ejected.bam",{"1":"chr1"},"helpme")
    elif args.filename == "gff":
        recontig.recontigGff(args.filename,"ejected.gff",{"1":"chr1"},"helpme")
    
    
if __name__ == "__main__":
    main()
