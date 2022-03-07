import pandas as pd
import glob
import recontig

mappings = {}

for x in glob.glob("fastas/*.fai"):
    org = x.split("/")[-1].split(".")[0]

    mappings[org] = {}
    mappings[org]["contigs"] = []

    mappings[org]["contigs"] = pd.read_table(x, names=["contigs","v2","v3","v4","v5"])["contigs"].to_list()

for from_org in list(mappings.keys()):
    
    for to_org in list(mappings.keys()):
        if from_org == to_org:
            continue
        
        mappings[from_org][to_org] = {}

        mapping = recontig.getContigMapping("recontig_mappings/{}_2_{}.mapping.txt".format(from_org, to_org))
        mappings[from_org][to_org]["recontig"] = mapping

        mapping = recontig.getContigMapping("re_mappings/{}_2_{}.mapping.txt".format(from_org, to_org))
        mappings[from_org][to_org]["regex"] = mapping

        mapping = recontig.getContigMapping("re_complex_mappings/{}_2_{}.mapping.txt".format(from_org, to_org))
        mappings[from_org][to_org]["regex_complex"] = mapping

        mapping = recontig.getContigMapping("checksum_mappings/{}_2_{}.mapping.txt".format(from_org, to_org))
        mappings[from_org][to_org]["raw_checksum"] = mapping

for from_org in list(mappings.keys()):
    print(from_org)
    for to_org in mappings[from_org].keys():
        print("\t" + to_org)

contig_to_row_mapping = {}

orgs = list(mappings.keys())
orgs.sort(reverse=True)

methods = ["recontig", "regex", "regex_complex", "raw_checksum"]

rows = []
i = 0

for from_org in orgs:
    for contig in mappings[from_org]["contigs"]:
        for method in methods:        
            already_found = False
            row_id = i
            row = {"method": method, from_org: 1}

            compatible_orgs = []
            contig_aliases = [contig]
            for to_org in orgs:
                if from_org == to_org:
                    continue
                
                # print(mappings[from_org])
                contigs = mappings[from_org][to_org][method]
                if contig in contigs:
                    compatible_orgs.append(to_org)
                    contig_aliases.append(contigs[contig])
            for org in compatible_orgs:
                row[org] = 1
            found = False
            for alias in contig_aliases:
                if alias + method in contig_to_row_mapping:
                    found = True
                    row_id = contig_to_row_mapping[alias + method]
                    break
            if found:
                rows[row_id].update(row)
                for alias in contig_aliases:
                    contig_to_row_mapping[alias + method] = row_id
            else:
                row["name"] = contig
                rows.append(row)
                for alias in contig_aliases:
                    contig_to_row_mapping[alias + method] = i
                i+=1
data = pd.DataFrame(rows)
data.fillna(0, inplace=True)
data = data[["name","method","UCSC","NCBI","gencode","ensembl"]]

print(data[(data["NCBI"] == 1) & (data["UCSC"] == 1)].shape)
data.to_csv("mapping_comparison.tsv", sep="\t", index=False)