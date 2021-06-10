module mapping;

import dhtslib.bgzf;
import std.array : split;

string[] BUILDS = [
    "BDGP6",
    "CanFam3",
    "GRCh37",
    "GRCh38",
    "GRCm37",
    "GRCm38",
    "GRCz10",
    "GRCz11",
    "JGI_4.2",
    "MEDAKA1",
    "R64-1-1",
    "Rnor_6.0",
    "WBcel235",
    "Xenopus_laevis_v2",
    "Xenopus_tropicalis_v9.1",
    "Zv9",
    "dm3",
    "galGal4",
    "galGal6",
    "rn5",
];

string[][] CONVERSIONS = [
    ["UCSC2ensembl", "ensembl2UCSC"], //BDGP6
    ["UCSC2ensembl"], //CanFam3
    [
        "NCBI2UCSC", "UCSC2ensembl", "UCSC2gencode", 
        "ensembl2UCSC", "ensembl2gencode", "gencode2UCSC", 
        "gencode2ensembl", 
    ], //GRCh37
    [
        "NCBI2ensembl", "RefSeq2UCSC", "UCSC2ensembl", 
        "UCSC2gencode", "ensembl2UCSC", "ensembl2gencode", 
        "gencode2UCSC", "gencode2ensembl", 
    ], //GRCh38
    [
        "UCSC2ensembl", "UCSC2gencode", "ensembl2UCSC", 
        "ensemblgencode", "gencode2UCSC", "gencode2ensembl",
        //^ not a typo
    ], //GRCm37
    [
        "UCSC2ensembl", "UCSC2gencode", "ensembl2UCSC", 
        "ensembl2gencode", "gencode2UCSC", "gencode2ensembl",
    ], //GRCm38
    ["UCSC2ensembl", "UCSC2gencode", "ensembl2UCSC", "gencode2UCSC",], //GRCz10
    ["UCSC2ensembl", "ensembl2UCSC"], //GRCz11
    ["UCSC2ensembl", "ensembl2UCSC"], //JGI_4.2
    ["UCSC2ensembl", "ensembl2UCSC"], //MEDAKA1
    ["UCSC2ensembl", "ensembl2UCSC"], //R64-1-1
    ["ensembl2UCSC"], //Rnor_6.0
    ["UCSC2ensembl", "ensembl2UCSC"], //WBcel235
    ["UCSC2xenbase", "xenbase2UCSC"], //Xenopus_laevis_v2
    ["UCSC2xenbase", "xenbase2UCSC"], //Xenopus_tropicalis_v9.1
    ["UCSC2ensembl", "ensembl2UCSC"], //Zv9
    ["UCSC2ensembl", "ensembl2UCSC"], //dm3
    ["UCSC2ensembl", "ensembl2UCSC"], //galGal4
    ["NCBI2UCSC", "UCSC2NCBI", "UCSC2ensembl", "ensembl2NCBI", "ensembl2UCSC"], //galGal6
    ["UCSC2ensembl", "ensembl2UCSC"], //rn5
];

auto getContigMappingFile(string build, string conversion)
{
    return BGZFile(
        "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/" ~
        build ~ "_" ~ conversion ~ ".txt"
        );
}

auto getContigMapping(string build, string conversion)
{
    string[string] mapping;
    auto file = getContigMappingFile(build, conversion);
    foreach (line; file.byLineCopy)
    {
        auto fields = line.split("\t");
        mapping[fields[0]] = fields[1];
    }
    return mapping;
}
