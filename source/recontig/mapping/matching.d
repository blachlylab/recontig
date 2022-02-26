module recontig.mapping.matching;
import std.digest.md;
import std.stdio;
import std.algorithm;
import std.range : iota;
import std.array : array;
import std.format : format;

import dhtslib.faidx;
import dhtslib.coordinates;
import htslib.hts_log;

import recontig.mapping.seq;
import recontig.mapping.checksum;

enum MatchByMasking{
    SM    = 1, // 001
    DN    = 2, // 010
    SM_DN = 3, // 011
    HM    = 4, // 100
    SM_HM = 5, // 101
    DN_HM = 6, // 110
    All   = 7, // 111
}

struct ContigMatcher
{
    IndexedFastaFile * fai1, fai2; 
    Checksum[string] fasta1Sums, fasta2Sums;
    string[][string] possiblyCompatible;
    string[string] compatible;
    string[] removed;

    this(IndexedFastaFile * fai1, IndexedFastaFile * fai2)
    {
        this.fai1 = fai1;
        this.fai2 = fai2;

        hts_log_info("recontig","Processing first fasta");
        this.fasta1Sums = this.processFasta(this.fai1);
        hts_log_info("recontig","Processing second fasta");
        this.fasta2Sums = this.processFasta(this.fai2);
    }

    /// get hashmap of Checksums for each contig in a fasta file
    Checksum[string] processFasta(IndexedFastaFile * fai){
        Checksum[string] fastaSums;

        /// loop over all contigs and calculate md5sum for fasta 1
        foreach (tid; iota(fai.nSeq))
        {
            auto chrom = fai.seqName(tid);
            /// get sections of fasta in chunks to keep memory usage low
            auto sums = checksumContig(fai, chrom);
            fastaSums[chrom] = sums;
        }

        return fastaSums;
    }
    
    /// Identify possible contig compatibility based on identical lengths
    void findCompatibleContigsByLength()
    {
        /// Find compatible contigs based on length
        foreach (contig1; fasta1Sums.byKey)
        {
            auto c1Len = fai1.seqLen(contig1);
            auto matched = false;
            foreach (contig2; fasta2Sums.byKey)
            {
                auto c2Len = fai2.seqLen(contig2);
                if(c1Len == c2Len){
                    possiblyCompatible[contig1] ~= contig2;
                    matched = true;
                }
            }
            if(!matched){
                hts_log_warning("recontig", "No possible mappings for contig %s based on length".format(contig1));
            }
        }
    }

    /// match together contigs that have identical raw checksums
    void matchContigsByRawCheckSum()
    {
        /// identify contigs with same raw checksum 
        foreach (contig1; possiblyCompatible.byKey.array)
        {
            auto matched = false;
            auto contig1sums = fasta1Sums[contig1];
            foreach (contig2; possiblyCompatible[contig1].dup)
            {
                auto contig2sums = fasta2Sums[contig2];
                if((contig1sums.hash == contig2sums.hash) && !removed.canFind(contig2)){
                    compatible[contig1] = contig2;
                    matched = true;
                    possiblyCompatible.remove(contig1);
                    removed ~= contig2;
                    break;
                }
            }
            if(!matched){
                hts_log_info("recontig", "No match found for contig %s using raw MD5 sum".format(contig1));
                hts_log_info("recontig", "We will attempt to modify the contig to find a match");
            }
        }
        this.dropRemoved;
    }

    /// remove any contigs from possiblyCompatible
    /// if they have been identified as compatible
    void dropRemoved()
    {
        /// remove matched contigs from hashmap
        foreach (contig1; possiblyCompatible.byKey.array)
        {
            foreach (contig2; possiblyCompatible[contig1].dup)
            {
                if(removed.canFind(contig2)){
                    auto idx = possiblyCompatible[contig1].countUntil(contig2);
                    possiblyCompatible[contig1] = possiblyCompatible[contig1].remove(idx);
                }
            }
        }
    }

    /// match together contigs that have identical checksums
    /// after modification using regions of:
    ///     degenerate nucleotides
    ///     soft-masking
    ///     hard-masking
    ///     or some combination of the above
    void matchContigsByMasking(MatchByMasking maskType)(){
        
        foreach (contig1; possiblyCompatible.byKey.array)
        {
            auto matched = false;
            foreach (i, contig2; possiblyCompatible[contig1].dup)
            {
                
                auto contig1NMaskedMd5 = this.applyMaskingToContig!maskType(contig1, contig2);
                auto contig2NMaskedMd5 = this.applyMaskingToContig!maskType(contig1, contig2, true);
                

                if((contig1NMaskedMd5 == contig2NMaskedMd5) && !removed.canFind(contig2)){
                    compatible[contig1] = contig2;
                    matched = true;
                    possiblyCompatible.remove(contig1);
                    removed ~= contig2;
                    break;
                }
            }
        }
        /// remove matched contigs from hashmap
        this.dropRemoved;
    }

    /// match together compatible contigs
    string[string] matchContigs(bool enforceMd5){

        this.findCompatibleContigsByLength();

        if(!enforceMd5)
        {
            hts_log_warning("recontig", "--no-enforce-md5 flag used");
            hts_log_warning("recontig", "Matching contigs purely on length");
            hts_log_warning("recontig", "mappings may not be the most accurate");
            /// remove matched contigs from hashmap
            foreach (contig1; possiblyCompatible.byKey)
            {                
                if(possiblyCompatible[contig1].length > 0)
                {
                    compatible[contig1] = possiblyCompatible[contig1][0];
                }
            }
            return compatible;
        }
        
        static foreach (i; 0..8)
        {
            static if(i == 0)
                this.matchContigsByRawCheckSum();
            else
                this.matchContigsByMasking!(cast(MatchByMasking) i)();            
        }

        /// report unmatched contigs
        foreach (contig1; possiblyCompatible.byKey)
        {
            hts_log_warning("recontig", "No match found for contig %s".format(contig1));
            if(possiblyCompatible[contig1].length > 0)
            {
                hts_log_warning("recontig", "Other unmatched contigs have the same length");
                hts_log_warning("recontig", "You may wish to use the --no-enforce-md5 flag");
            }
        }
        
        return compatible;
    }

    /// collect degenerate regions for both contigs
    auto degenerateRegions(string chrom1, string chrom2)
    {
        return (
            this.fasta1Sums[chrom1].degenerateRegions ~ 
            this.fasta2Sums[chrom2].degenerateRegions
        ).sort!"a.start < b.start"
        .chunkBy!"a.isOverlap(b)"
        .map!(x => x.fold!"a | b")
        .array;
    }
    
    /// collect soft-masked regions for both contigs
    auto softMaskedRegions(string chrom1, string chrom2)
    {
        return (
            this.fasta1Sums[chrom1].softMaskedRegions ~ 
            this.fasta2Sums[chrom2].softMaskedRegions
        ).sort!"a.start < b.start"
        .chunkBy!"a.isOverlap(b)"
        .map!(x => x.fold!"a | b")
        .array;
    }

    /// collect hard-masked regions for both contigs
    auto hardMaskedRegions(string chrom1, string chrom2)
    {
        return (
            this.fasta1Sums[chrom1].hardMaskedRegions ~ 
            this.fasta2Sums[chrom2].hardMaskedRegions
        ).sort!"a.start < b.start"
        .chunkBy!"a.isOverlap(b)"
        .map!(x => x.fold!"a | b")
        .array;
    }

    /// modify contigs using regions to modify MD5 sum
    /// Allows us to match together contigs with differential:
    /// degenerate nucleotides
    /// soft-masking
    /// hard-masking
    auto applyMaskingToContig(MatchByMasking maskType)(string chrom1, string chrom2, bool reverse = false)
    {    
        MD5 md5sum;

        assert(fai1.seqLen(chrom1) == fai2.seqLen(chrom2));

        md5sum.start();
        /// get sections of fasta in chunks to keep memory usage low
        static if(maskType & 1)
            auto smUnion = this.softMaskedRegions(chrom1, chrom2);
        static if(maskType & 2)
            auto dgUnion = this.degenerateRegions(chrom1, chrom2);
        static if(maskType & 4)
            auto hmUnion = this.hardMaskedRegions(chrom1, chrom2);

        foreach (i; iota(0, fai1.seqLen(chrom1), 4_000_000))
        {
            auto end = i + 4_000_000 > fai1.seqLen(chrom1) ? fai1.seqLen(chrom1) : 4_000_000;
            auto coords = ZBHO(i, i + end);
            auto seq1 = reverse ? (*fai2)[chrom2, coords].dup : (*fai1)[chrom1, coords].dup;
            
            static if(maskType & 1)
                convertSoftMaskedRegions(seq1, smUnion, coords);
            static if(maskType & 2){
                assert(fai2);
                auto seq2 = reverse ? (*fai1)[chrom1, coords].dup : (*fai2)[chrom2, coords].dup;
                convertDegenerateRegions(seq1, seq2, dgUnion, coords);
            }static if(maskType & 4)
                hardMaskRegions(seq1, hmUnion, coords);

            md5sum.put(cast(const(ubyte)[]) seq1);
        }
        return toHexString(md5sum.finish()).idup;
    }

}

unittest
{
    auto fai1 = IndexedFastaFile("tests/data/baseline.fa");
    auto fai2 = IndexedFastaFile("tests/data/test1.fa");

    ContigMatcher cm = ContigMatcher(&fai1, &fai2);
    
    auto cs1 = cm.applyMaskingToContig!(MatchByMasking.DN)("chrA", "A_A");
    auto cs2 = cm.applyMaskingToContig!(MatchByMasking.DN)("chrA", "A_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.HM)("chrB", "B_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.HM)("chrB", "B_A", true);
    assert(cs1 == cs2);
    
    cs1 = cm.applyMaskingToContig!(MatchByMasking.DN)("chrC", "C_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.DN)("chrC", "C_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.SM)("chrD", "D_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.SM)("chrD", "D_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.DN_HM)("chrBC", "BC_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.DN_HM)("chrBC", "BC_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.SM_HM)("chrBD", "BD_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.SM_HM)("chrBD", "BD_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.SM_DN)("chrCD", "CD_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.SM_DN)("chrCD", "CD_A", true);
    assert(cs1 == cs2);

    cs1 = cm.applyMaskingToContig!(MatchByMasking.All)("chrBCD", "BCD_A");
    cs2 = cm.applyMaskingToContig!(MatchByMasking.All)("chrBCD", "BCD_A", true);
    assert(cs1 == cs2);
}

unittest
{
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    auto fai1 = IndexedFastaFile("tests/data/baseline.fa");
    auto fai2 = IndexedFastaFile("tests/data/test1.fa");

    ContigMatcher cm = ContigMatcher(&fai1, &fai2);

    auto mapping = cm.matchContigs(true);
    assert(mapping["chrA"]     == "A_A");
    assert(mapping["chrB"]     == "B_A");
    assert(mapping["chrC"]     == "C_A");
    assert(mapping["chrD"]     == "D_A");
    assert(mapping["chrBC"]   == "BC_A");
    assert(mapping["chrBD"]   == "BD_A");
    assert(mapping["chrBCD"] == "BCD_A");

    fai2 = IndexedFastaFile("tests/data/test2.fa");
    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_B");
    assert(mapping["chrB"]     == "B_B");
    assert(mapping["chrC"]     == "C_B");
    assert(mapping["chrD"]     == "D_B");
    assert(mapping["chrBC"]   == "BC_B");
    assert(mapping["chrBD"]   == "BD_B");
    assert(mapping["chrBCD"] == "BCD_B");

    fai2 = IndexedFastaFile("tests/data/test3.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_C");
    assert(mapping["chrB"]     == "B_C");
    assert(mapping["chrC"]     == "C_C");
    assert(mapping["chrD"]     == "D_C");
    assert(mapping["chrBC"]   == "BC_C");
    assert(mapping["chrBD"]   == "BD_C");
    assert(mapping["chrBCD"] == "BCD_C");

    fai2 = IndexedFastaFile("tests/data/test4.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_D");
    assert(mapping["chrB"]     == "B_D");
    assert(mapping["chrC"]     == "C_D");
    assert(mapping["chrD"]     == "D_D");
    assert(mapping["chrBC"]   == "BC_D");
    assert(mapping["chrBD"]   == "BD_D");
    assert(mapping["chrBCD"] == "BCD_D");

    fai2 = IndexedFastaFile("tests/data/test5.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);

    assert(mapping["chrA"]     == "A_BC");
    assert(mapping["chrB"]     == "B_BC");
    assert(mapping["chrC"]     == "C_BC");
    assert(mapping["chrD"]     == "D_BC");
    assert(mapping["chrBC"]   == "BC_BC");
    assert(mapping["chrBD"]   == "BD_BC");
    assert(mapping["chrBCD"] == "BCD_BC");

    fai2 = IndexedFastaFile("tests/data/test6.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_BD");
    assert(mapping["chrB"]     == "B_BD");
    assert(mapping["chrC"]     == "C_BD");
    assert(mapping["chrD"]     == "D_BD");
    assert(mapping["chrBC"]   == "BC_BD");
    assert(mapping["chrBD"]   == "BD_BD");
    assert(mapping["chrBCD"] == "BCD_BD");

    fai2 = IndexedFastaFile("tests/data/test7.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_CD");
    assert(mapping["chrB"]     == "B_CD");
    assert(mapping["chrC"]     == "C_CD");
    assert(mapping["chrD"]     == "D_CD");
    assert(mapping["chrBC"]   == "BC_CD");
    assert(mapping["chrBD"]   == "BD_CD");
    assert(mapping["chrBCD"] == "BCD_CD");

    fai2 = IndexedFastaFile("tests/data/test8.fa");

    cm = ContigMatcher(&fai1, &fai2);

    mapping = cm.matchContigs(true);
    
    assert(mapping["chrA"]     == "A_BCD");
    assert(mapping["chrB"]     == "B_BCD");
    assert(mapping["chrC"]     == "C_BCD");
    assert(mapping["chrD"]     == "D_BCD");
    assert(mapping["chrBC"]   == "BC_BCD");
    assert(mapping["chrBD"]   == "BD_BCD");
    assert(mapping["chrBCD"] == "BCD_BCD");
}