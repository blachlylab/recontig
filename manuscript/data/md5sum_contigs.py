import hashlib
import argparse

def md5_contigs(f):
    contig_sums = {}
    md5 = hashlib.md5()
    contig = ""
    # i = 0
    with open(f, "r", encoding='utf-8') as a_file:
        for line in a_file:
            line = line.strip()
            if contig == "":
                contig = line[1::].split(" ")[0]
                # print(contig)
                continue
            
            if line.startswith(">"):
                # i = 0
                contig_sums[md5.hexdigest()] = contig
                md5 = hashlib.md5()
                contig = line[1::].split(" ")[0]
                # print(contig)
                continue
                  
            # if i == 10:
            #     continue
            # print(line+"end")
            md5.update(line.encode('utf-8'))
            # i+=1
        contig_sums[md5.hexdigest()] = contig
    return contig_sums
        
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta1")
    parser.add_argument("fasta2")
    args = parser.parse_args()

    sums1 = md5_contigs(args.fasta1)
    sums2 = md5_contigs(args.fasta2)

    for k in sums1:
        if k in sums2:
            print("{}\t{}".format(sums1[k],sums2[k]))
        

