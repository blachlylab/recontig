module recontig.mapping;

import std.stdio;
import std.digest.md;
import std.array : split;
import std.range : iota;
import std.algorithm : map, sort;
import std.array : array, split;
import std.format : format;

import dhtslib.bgzf;
import dhtslib.faidx;
import dhtslib.coordinates;
import htslib.hts : seq_nt16_table, seq_nt16_int;
import htslib.hts_log;

/// Builds specified on dpryan79's ChromosomeMappings github
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

/// conversions for each build specified on dpryan79's ChromosomeMappings github
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


/// download a mapping from dpryan79's ChromosomeMappings github
auto getDpryan79ContigMapping(string build, string conversion)
{
    return BGZFile(
        "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/" ~
        build ~ "_" ~ conversion ~ ".txt"
        ).convertMappingToHashMap;
    
}

/// load a contig mapping file
auto getContigMapping(string fn)
{
    return BGZFile(fn).convertMappingToHashMap;
}

/// load a contig mapping file into a hashmap
private auto convertMappingToHashMap(BGZFile file)
{
    string[string] mapping;
    foreach (line; file.byLineCopy)
    {
        auto fields = line.split("\t");
        if(fields[0] == "" || fields[1] == "")
            hts_log_info("recontig", "line %s in mapping file skipped".format(line));
        mapping[fields[0]] = fields[1];
    }
    return mapping;
}

/// make a contig mapping file from two faidx'd fasta files
/// Only chromosomes/contigs that generate the 
/// same md5sum are considered equal and 
/// output as a mapping
/// fasta sequences are corrected by
/// making all nucleotides uppercase and 
/// converting all degenerate nucleotides to N
void makeMapping(string fa1, string fa2)
{

    // load faidx'd fasta files
    auto fai1 = IndexedFastaFile(fa1);
    auto fai2 = IndexedFastaFile(fa2);

    // create hashmaps
    string[string] fasta1Sums;
    string[string] fasta2Sums;
    stderr.writefln("[recontig]: Processing first fasta");

    /// loop over all contigs and calculate md5sum for fasta 1
    foreach (tid; iota(fai1.nSeq))
    {
        MD5 md5sum;
        md5sum.start();
        auto chrom = fai1.seqName(tid);
        /// get sections of fasta in chunks to keep memory usage low
        foreach (i; iota(0, fai1.seqLen(chrom), 100_000))
        {
            auto end = i + 100_000 > fai1.seqLen(chrom) ? fai1.seqLen(chrom) : 100_000;
            auto coords = ZBHO(i, i + end);

            auto bases = fai1[chrom, coords].cleanFasta;

            md5sum.put(cast(const(ubyte)[])bases);
        }
        fasta1Sums[toHexString(md5sum.finish()).idup] = chrom;
    }
    stderr.writefln("[recontig]: Processing second fasta");
    /// loop over all contigs and calculate md5sum for fasta 2
    foreach (tid; iota(fai2.nSeq))
    {
        MD5 md5sum;
        md5sum.start();
        auto chrom = fai2.seqName(tid);
        /// get sections of fasta in chunks to keep memory usage low
        foreach (i; iota(0, fai2.seqLen(chrom), 100_000))
        {
            auto end = i + 100_000 > fai2.seqLen(chrom) ? fai2.seqLen(chrom) : 100_000;
            auto coords = ZBHO(i, i + end);
            // auto region1 = ChromCoordinates(chrom1, coords);
            // auto region2 = ChromCoordinates(chrom2, coords);
            auto bases = fai2[chrom, coords].cleanFasta;
            md5sum.put(cast(const(ubyte)[])bases);
        }
        fasta2Sums[toHexString(md5sum.finish()).idup] = chrom;
        
    }
    /// intersects md5sums between the two fasta's to create a contig mapping
    foreach (item; fasta1Sums.byKeyValue.array.sort!((a, b) => a.value < b.value))
    {
        if(item.key in fasta2Sums){
            writefln("%s\t%s",fasta1Sums[item.key],fasta2Sums[item.key]);
        }
    }    
}

/// used along with seq_nt16_int & seq_nt16_table to convert
/// all lowercase nucleotides and non ACGTN nucleotides
/// to N
char[5] iupacToACGTN = ['A','C','G','T','N'];

/// Convert nucleotide sequences to upper case and
/// remove IUPAC/degenerate nucleotides
pragma(inline, true) auto cleanFasta(string seq)
{
    return seq.map!(x => iupacToACGTN[seq_nt16_int[seq_nt16_table[x]]]).array;
}