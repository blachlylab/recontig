module recontig.mapping.generate;

import std.digest.md;
import std.stdio;
import std.algorithm;
import std.range : iota, zip;
import std.array : array;
import std.format : format;

import dhtslib.faidx;
import dhtslib.coordinates;
import htslib.hts_log;

import recontig.mapping.seq;
import recontig.mapping.checksum;
import recontig.mapping.matching;

/// make a contig mapping file from two faidx'd fasta files
/// output to specified file
void makeMapping(string fa1, string fa2, string fo = "-", bool enforceMd5 = true)
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

/// make a contig mapping file from two faidx'd fasta files
string[string] makeMapping(string fa1, string fa2, bool enforceMd5 = false)
{
    string[string] ret; 
    // load faidx'd fasta files
    auto fai1 = IndexedFastaFile(fa1);
    auto fai2 = IndexedFastaFile(fa2);

    // create hashmaps
    ContigMatcher cm = ContigMatcher(&fai1, &fai2);
    
    /// intersects md5sums between the two fasta's to create a contig mapping
    return cm.matchContigs(enforceMd5);
}

