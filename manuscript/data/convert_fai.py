import re
import argparse
import pandas as pd


# can convert major contigs from human reference builds for
# gencode and UCSC to ensembl format, basically removes chr prefix
# 
# i.e chr1 -> 1 
# 
# can easily be replicated with sed as
# sed -E 's/^chr([0-9XYMUn]+)/\1/g'
regex_naive_conversion = ["^chr([0-9XYMUn]+)","\g<1>"]

# can convert major contigs and alternate contigs from human reference builds for
# gencode, UCSC, and NCBI to ensembl format
# 
# harder to replicate into sed
regex_complex_conversions = [["^chr([0-9XYMUn]+)","\g<1>"],["^[0-9XYM]+_([0-9a-zA-Z]+)v([0-9]).*","\g<1>.\g<2>"],["^M$","MT"],["^NC_[0]+([0-9]+)\.[0-9]+","\g<1>"],["^23$","Y"],["^24$","X"],["^23$","Y"],["^12920$","M"]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--complex", action="store_true")
    parser.add_argument("fai1")
    parser.add_argument("fai2")
    args = parser.parse_args()

    contigs1 = pd.read_table(args.fai1, names=["contigs","v2","v3","v4","v5"])["contigs"].to_list()
    contigs2 = pd.read_table(args.fai2, names=["contigs","v2","v3","v4","v5"])["contigs"].to_list()

    fai1 = {}
    fai2 = {}

    converted1 = {}
    converted2 = {}

    for i, contig in enumerate(contigs1):
        fai1[i] = contig
        if args.complex:
            for reg in regex_complex_conversions:
                contig = re.sub(reg[0], reg[1], contig)
        else:
            contig = re.sub(regex_naive_conversion[0], regex_naive_conversion[1], contig)
        converted1[contig] = i

    for i, contig in enumerate(contigs2):
        fai2[i] = contig
        if args.complex:
            for reg in regex_complex_conversions:
                contig = re.sub(reg[0], reg[1], contig)
        else:
            contig = re.sub(regex_naive_conversion[0], regex_naive_conversion[1], contig)
        converted2[contig] = i
        
    for k in converted1:
        if k in converted2:
            print("{}\t{}".format(fai1[converted1[k]], fai2[converted2[k]]))


        


