module recontig.mapping;

import std.stdio;
import std.digest.md;
import std.array : split;
import std.range : iota;
import std.algorithm : map, sort;
import std.array : array, split;
import std.format : format;
import std.conv : to;

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
    if(build =="" || conversion == ""){
        throw new Exception("Error: if not using a mapping file you must provide a valid build and conversion.");
    }
    // validate build
    bool buildFound;
    ulong buildIdx = 0;
    foreach(i, b;BUILDS){
        if(b == build) buildFound = true, buildIdx = i;
    }
    if(!buildFound){
        throw new Exception("Error: Please use a valid build: " ~ BUILDS.to!string);
    }

    //validate conversion
    bool convFound;
    ulong convIdx = 0;
    foreach(i, c;CONVERSIONS[buildIdx]){
        if(c == conversion) convFound = true, convIdx = i;
    }
    if(!convFound){
        throw new Exception("Error: Please use a valid conversion: " ~ CONVERSIONS[buildIdx].to!string);
    }
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
void makeMapping(string fa1, string fa2, string fo, bool enforceMd5 = false)
{
    File f;
    if(fo == "" || fo == "-"){
        f = stdout;
    }else{
        f = File(fo, "w");
    }
    auto mapping = makeMapping(fa1, fa2, enforceMd5);
    /// intersects md5sums between the two fasta's to create a contig mapping
    foreach (item; mapping.byKeyValue.array.sort!((a, b) => a.key < b.key))
    {
        f.writefln("%s\t%s",item.key,item.value);
    }    
}

string[string] makeMapping(string fa1, string fa2, bool enforceMd5 = false)
{
    string[string] ret; 
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
    foreach (contig1; fasta1Sums.byKeyValue.array.sort!((a, b) => a.value < b.value))
    {
        auto c1Len = fai1.seqLen(contig1.value);
        auto matched = false;
        foreach (contig2; fasta2Sums.byKeyValue.array.sort!((a, b) => a.value < b.value))
        {
            auto c2Len = fai2.seqLen(contig2.value);
            if(contig1.key == contig2.key){
                assert(c1Len == c2Len);
                ret[fasta1Sums[contig1.key]] = fasta2Sums[contig1.key];
                matched = true;
                break;
            }else if(c1Len == c2Len){
                if(!enforceMd5){
                    hts_log_warning("recontig", "contigs %s and %s lengths match but their MD5 checksums do not.".format(contig1.value, contig2.value));
                    ret[fasta1Sums[contig1.key]] = fasta2Sums[contig2.key];
                    matched = true;
                    break;
                }else{
                    hts_log_warning("recontig", "contigs %s and %s lengths match but their MD5 checksums do not.".format(contig1, contig2));
                    hts_log_warning("recontig", "--enforce-md5sum flag was used. This pair will not be output");
                }   
            }
        }
        if(!matched){
            hts_log_warning("recontig", "No match found for contig %s".format(contig1.value));
        }
    }
    return ret;
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