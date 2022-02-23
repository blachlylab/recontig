module recontig.mapping.checksum;

import std.digest.md;
import std.range: iota;

import dhtslib.coordinates;
import dhtslib.faidx;

import recontig.mapping.seq;

struct ChecksumBuilder
{

    MD5 hash;

    ZBHO[] hardMaskedRegions;
    bool lastHMRegionFinished = true;

    ZBHO[] softMaskedRegions;
    bool lastSMRegionFinished = true;

    ZBHO[] degenerateRegions;
    bool lastDRegionFinished = true;

    /// initialize MD5 sum
    void initialize()
    {
        this.hash.start();
    }

    /// append new seq data
    /// also record relevent regions
    void append(string seq, ZBHO coords)
    {
        hash.put(cast(const(ubyte)[]) seq);

        getHardMaskedRegions(seq, coords.start, this.hardMaskedRegions, &this.lastHMRegionFinished);
        getSoftMaskedRegions(seq, coords.start, this.softMaskedRegions, &this.lastSMRegionFinished);
        getDegerateRegions(seq, coords.start, this.degenerateRegions, &this.lastDRegionFinished);
    }

    /// finalize md5 sum
    string finalize()
    {
        string ret = toHexString(this.hash.finish()).idup;
        return ret;
    }
}

struct Checksum
{
    string hash;      

    ZBHO[] hardMaskedRegions;
    ZBHO[] softMaskedRegions;
    ZBHO[] degenerateRegions;

    this(ChecksumBuilder cs)
    {
        this.hash = cs.finalize;
        
        this.hardMaskedRegions    = cs.hardMaskedRegions;
        this.softMaskedRegions    = cs.softMaskedRegions;
        this.degenerateRegions    = cs.degenerateRegions;
    }

    /// allow for use in a hashmap
    size_t toHash() const @nogc @safe pure nothrow
    {
        return hash.hashOf();
    }
}

/// Build checksum for contig
auto checksumContig(IndexedFastaFile * fai, string chrom)
{
    ChecksumBuilder sums;
    sums.initialize;

    /// get sections of fasta in chunks to keep memory usage low
    foreach (i; iota(0, fai.seqLen(chrom), 100_000))
    {
        auto end = i + 100_000 > fai.seqLen(chrom) ? fai.seqLen(chrom) : 100_000;
        auto coords = ZBHO(i, i + end);
        auto seq = (*fai)[chrom, coords];

        sums.append(seq, coords);
    }
    
    Checksum cs = Checksum(sums);

    return cs;
}

unittest
{
    import std.range: repeat, takeExactly;
    import std.array: array;
    import std.stdio;

    auto fai = IndexedFastaFile("tests/data/baseline.fa");
    
    auto cs = checksumContig(&fai, "chrA");
    assert(cs.hash == "CE80DFA04030E641A9B4CD9B7C4E33D2");
    assert(cs.hardMaskedRegions == []);
    assert(cs.softMaskedRegions == []);
    assert(cs.degenerateRegions == []);


    cs = checksumContig(&fai, "chrB");
    assert(cs.hash == "5ED8C429355FA340719CCBC0AD6B8F39");
    assert(cs.hardMaskedRegions == [ZBHO(600,720)]);
    assert(cs.softMaskedRegions == []);
    assert(cs.degenerateRegions == []);
    foreach (ZBHO key; cs.hardMaskedRegions){
        assert(fai["chrB", key] == 'N'.repeat.takeExactly(key.size).array.idup);
    }

    cs = checksumContig(&fai, "chrC");
    assert(cs.hash == "F235D6D02BB91D79DD0AFC4DEF8C1B2B");
    assert(cs.hardMaskedRegions == []);
    assert(cs.softMaskedRegions == []);
    assert(cs.degenerateRegions == [ZBHO(480,540)]);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrC", key].convertAllToHardMask == 'N'.repeat.takeExactly(key.size).array.idup);

    cs = checksumContig(&fai, "chrD");
    assert(cs.hash == "5E8531C95375E4992821B8DB126664EA");
    assert(cs.hardMaskedRegions == []);
    assert(cs.softMaskedRegions == [ZBHO(480,540)]);
    assert(cs.degenerateRegions == []);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrD", key].convertAllToHardMask == 'N'.repeat.takeExactly(key.size).array.idup);

    cs = checksumContig(&fai, "chrBC");
    assert(cs.hash == "7AB1E2A25BC23D1384453197BE103452");
    assert(cs.hardMaskedRegions == [ZBHO(0,29), ZBHO(928,960)]);
    assert(cs.softMaskedRegions == []);
    assert(cs.degenerateRegions == [ZBHO(480,540)]);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrBC", key].convertAllToHardMask == 'N'.repeat.takeExactly(key.size).array.idup);

    cs = checksumContig(&fai, "chrBD");
    assert(cs.hash == "212FA1B70E1FAAC451D35A2711315006");
    assert(cs.hardMaskedRegions == [ZBHO(0, 29), ZBHO(808, 840)]);
    assert(cs.softMaskedRegions == [ZBHO(360,420)]);
    assert(cs.degenerateRegions == []);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrBD", key].convertAllToHardMask == 'N'.repeat.takeExactly(key.size).array.idup);

    cs = checksumContig(&fai, "chrCD");
    assert(cs.hash == "97E9C4D52A9F237B59CCA81454D035F0");
    assert(cs.hardMaskedRegions == []);
    assert(cs.softMaskedRegions == [ZBHO(360,420)]);
    assert(cs.degenerateRegions == [ZBHO(480,540)]);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrCD", key].convertAllToHardMask== 'N'.repeat.takeExactly(key.size).array.idup);
    
    cs = checksumContig(&fai, "chrBCD");
    assert(cs.hash == "A89584168A19218755691A079F69AB8D");
    assert(cs.hardMaskedRegions == [ZBHO(0, 29), ZBHO(748, 780)]);
    assert(cs.softMaskedRegions == [ZBHO(300, 360)]);
    assert(cs.degenerateRegions == [ZBHO(420, 480)]);
    foreach (ZBHO key; cs.hardMaskedRegions)
        assert(fai["chrBCD", key].convertAllToHardMask == 'N'.repeat.takeExactly(key.size).array.idup);
}
