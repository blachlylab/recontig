import pandas as pd
import glob
from pkg_resources import compatible_platforms
import recontig

mappings = {}

for x in glob.glob("fastas/*.fai"):
    org = x.split("/")[-1].split(".")[0]

    mappings[org] = {}
    mappings[org]["contigs"] = set()

    mappings[org]["contigs"].update(pd.read_table(x, names=["contigs","v2","v3","v4","v5"])["contigs"].to_list())

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




orgs = list(mappings.keys())
orgs.sort(reverse=True)

methods = ["recontig", "regex", "regex_complex", "raw_checksum"]

unique_contig_data = {}
row_to_contigs = {}

unique_contig = set()
for org in orgs:
    unique_contig.update(mappings[org]["contigs"])


rows = []
num_uniq_contigs = 0
already_used = {}

for contig in unique_contig:
    compatible_contigs = set()
    compatible_contigs.add(contig)

    # collect all contigs that are compatible
    # with this contig entry
    while True:
        contig_list = list(compatible_contigs)
        prior_len = len(contig_list)
        for c in contig_list:
            for method in methods:
                for from_org in orgs:
                    for to_org in orgs:
                        if from_org == to_org:
                            continue
                        other = mappings[from_org][to_org][method]
                        if c in other:
                            if other[c] not in compatible_contigs:
                                compatible_contigs.add(other[c])

        if prior_len == len(compatible_contigs):
            break

    found = False
    contig_id = num_uniq_contigs
    # if this contig already exists with other names 
    # use that entry
    for c in compatible_contigs:
        if c in already_used:
            contig_id = already_used[c]
            found = True
            break
    
    # if this contig doesn't already exists with other names 
    # create entry
    if not found:
        unique_contig_data[contig_id] = {}
        for method in methods:
            
            unique_contig_data[contig_id][method] = {}
            unique_contig_data[contig_id][method] = {"contigs": set()}
        num_uniq_contigs+=1

    for method in methods:
        unique_contig_data[contig_id][method]["contigs"].update(compatible_contigs)

    # mark as already exists
    for x in compatible_contigs:
        already_used[c] = contig_id

# loop over compatible contig entrys
for id,contigs in unique_contig_data.items():
    # loop over methods
    for method,contig in contigs.items():
        contig["orgs"] = set()
        # for each entry collect org names
        for c in contig["contigs"]:    
            for from_org in orgs:
                for to_org in orgs:
                    if from_org == to_org:
                        continue
                    if c in mappings[from_org]["contigs"]:
                        other = mappings[from_org][to_org][method]
                        if c in other:
                            contig["orgs"].add(from_org)
                            contig["orgs"].add(to_org)
    # convert sets to lists
    for method,contig in contigs.items():
        contig["contigs"] = list(contig["contigs"])
        contig["orgs"] = list(contig["orgs"])

# create matrix from dictionary for pandas
for id,contigs in unique_contig_data.items():
    for method,v in contigs.items():
        data = {"name": v["contigs"][0], "method": method}
        for org in orgs:
            data[org] = 0
        for org in v["orgs"]:
            data[org] = 1
        rows.append(data)

data = pd.DataFrame(rows)
data.fillna(0, inplace=True)
data = data[["name","method","UCSC","NCBI","gencode","ensembl"]]

data.to_csv("mapping_comparison.tsv", sep="\t", index=False)

# loop over compatible contig entrys
for id,contigs in unique_contig_data.items():
    # loop over methods
    for method,contig in contigs.items():
        if len(contig["orgs"]) == 0:
            for org in orgs:
                if contig["contigs"][0] in mappings[org]["contigs"]:
                    contig["orgs"].append(org)
                    break
rows = []
# create matrix from dictionary for pandas
for id,contigs in unique_contig_data.items():
    for method,v in contigs.items():
        data = {"name": v["contigs"][0], "method": method}
        for org in orgs:
            data[org] = 0
        for org in v["orgs"]:
            data[org] = 1
        rows.append(data)

data = pd.DataFrame(rows)
data.fillna(0, inplace=True)
data = data[["name","method","UCSC","NCBI","gencode","ensembl"]]
data["sum"] = data[["UCSC","NCBI","gencode","ensembl"]].sum(axis=1)
data.sort_values("sum", inplace=True, ascending=True)
data.drop("sum", axis=1, inplace=True)


data.to_csv("mapping_comparison_with_singletons.tsv", sep="\t", index=False)

counts = {}
for from_org in orgs:
    counts[from_org] = {}
    for to_org in orgs:
        if from_org == to_org:
            continue
        mappings[from_org]["contigs"]
        a_and_b = len(mappings[from_org][to_org]["recontig"])
        a_diff_b = len(mappings[from_org]["contigs"] - set(mappings[from_org][to_org]["recontig"].keys()))
        b_diff_a = len(mappings[to_org]["contigs"] - set(mappings[from_org][to_org]["recontig"].values()))
        counts[from_org][to_org] = (a_and_b, a_diff_b, b_diff_a)
data = pd.DataFrame(counts).sort_index()
data = data[data.index.tolist()]
data.fillna("", inplace=True)

import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(10,3))
ax=plt.subplot(111)
ax.axis('off')
c = data.shape[1]
r = data.shape[0]
import math
def print_tuple(x):
    ret = "( "
    start = 1
    if type(x) == tuple:
        for i in x:
            num = str(i)
            digits = 0
            if i > 0:
                digits = int(math.log10(i))+1
            elif i == 0:
                digits = 1
            if digits != 3:
                ret +=(" "*((3-(digits%3))+2))
            ret += str(i)
            ret+=", "
            start = len(ret)
        ret = ret[:len(ret)-2]+" "
        ret += ")"
        return ret
    else:
        return x

values = data.applymap(print_tuple)
table_data = np.hstack([[[x] for x in [""] + data.index.tolist()], np.vstack([data.columns, values])])
cellColours = [['none'] + ['lightgray']*c] + [['lightgray'] + ['none']*c]*r
ax.table(cellText=table_data, cellColours=cellColours,  bbox=[0,0,1,1])
plt.savefig("table1.pdf")

import matplotlib as mpl

def mergecells(table, cells):
    '''
    Merge N matplotlib.Table cells

    Parameters
    -----------
    table: matplotlib.Table
        the table
    cells: list[set]
        list of sets od the table coordinates
        - example: [(0,1), (0,0), (0,2)]

    Notes
    ------
    https://stackoverflow.com/a/53819765/12684122
    '''
    cells_array = [np.asarray(c) for c in cells]
    h = np.array([cells_array[i+1][0] - cells_array[i][0] for i in range(len(cells_array) - 1)])
    v = np.array([cells_array[i+1][1] - cells_array[i][1] for i in range(len(cells_array) - 1)])

    # if it's a horizontal merge, all values for `h` are 0
    if not np.any(h):
        # sort by horizontal coord
        cells = np.array(sorted(list(cells), key=lambda v: v[1]))
        edges = ['BTL'] + ['BT' for i in range(len(cells) - 2)] + ['BTR']
    elif not np.any(v):
        cells = np.array(sorted(list(cells), key=lambda h: h[0]))
        edges = ['TRL'] + ['RL' for i in range(len(cells) - 2)] + ['BRL']
    else:
        raise ValueError("Only horizontal and vertical merges allowed")

    for cell, e in zip(cells, edges):
        table[cell[0], cell[1]].visible_edges = e
        
    txts = [table[cell[0], cell[1]].get_text() for cell in cells]
    tpos = [np.array(t.get_position()) for t in txts]

    # transpose the text of the left cell
    trans = (tpos[-1] - tpos[0])/2
    # didn't had to check for ha because I only want ha='center'
    txts[0].set_transform(mpl.transforms.Affine2D().translate(*trans))
    for txt in txts[1:]:
        txt.set_visible(False)

counts = {}
counts["Source ∩ Destination"] = {}
counts["Source - Destination"] = {}
counts["Destination - Source"] = {}
for from_org in orgs:
    counts["Source ∩ Destination"][from_org] = {}
    counts["Source - Destination"][from_org] = {}
    counts["Destination - Source"][from_org] = {}
    for to_org in orgs:
        if from_org == to_org:
            continue
        mappings[from_org]["contigs"]
        a_and_b = len(mappings[from_org][to_org]["recontig"])
        a_diff_b = len(mappings[from_org]["contigs"] - set(mappings[from_org][to_org]["recontig"].keys()))
        b_diff_a = len(mappings[to_org]["contigs"] - set(mappings[from_org][to_org]["recontig"].values()))
        counts["Source ∩ Destination"][from_org][to_org] = a_and_b
        counts["Source - Destination"][from_org][to_org] = a_diff_b
        counts["Destination - Source"][from_org][to_org] = b_diff_a

reform = {(outerKey, innerKey): values for outerKey, innerDict in counts.items() for innerKey, values in innerDict.items()}

data = pd.DataFrame(reform)
data = data[data.sort_index(level=[0,1], ascending=[False,True], axis=1, ).columns]
data = data.sort_index()

data.fillna("", inplace=True)
fig = plt.figure(figsize=(8,3))
ax=fig.gca()
ax.axis('off')
c = data.shape[1] + 1
r = data.shape[0]

toplevel = data.columns.get_level_values(0).tolist()
top_row = ['',toplevel[0],'','','',toplevel[4],'','','',toplevel[8],'','','']
next_row = [''] + data.columns.get_level_values(1).tolist()
values_and_index = np.hstack([[[x] for x in data.index.tolist()], data.values])
table_data = np.vstack([[top_row, next_row], values_and_index])

ax.table(cellColours=[['lightgray']] + [['none']] + [['none']], bbox=[0,0,1,1])
# print(cellColours)
table = ax.table(cellText=table_data, cellColours=[['none']*c]*(2 + r),  bbox=[0,0,1,1])

fig.canvas.draw()
# mergecells(table, (1,0), (0,0))

mergecells(table, [(0,1), (0,2), (0,3), (0,4)])
mergecells(table, [(0,5), (0,6), (0,7), (0,8)])
mergecells(table, [(0,9), (0,10), (0,11), (0,12)])
table.auto_set_font_size(False)
table.set_fontsize(6)

plt.show()