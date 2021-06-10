module recontig.bam;

import std.utf : toUTFz;
import std.conv : to;

import dhtslib.sam;
import htslib.sam : sam_hdr_remove_line_id;

/// recontig bam/sam to bam file
alias recontigBam = recontigBamImpl!true;

/// recontig bam/sam to sam file
alias recontigSam = recontigBamImpl!false; 


/// recontig BAM/SAM file
/// makes edits to header SQ records
/// changes bam1_t tid and mate tid to reflect new contig names
void recontigBamImpl(bool outputBam)(string fn, string ejectedfn, string[string] mapping, string argStr)
{
    // Open new bam readers
    auto bamr = SAMReader(fn);
    auto ejectedbamw = SAMWriter(ejectedfn, bamr.header, SAMWriterTypes.SAM);

    // Header manipulation
    auto newHeader = bamr.header.dup;

    // remove old SQ records
    foreach (key; mapping.byKey)
    {
        sam_hdr_remove_line_id(newHeader.h, toUTFz!(char *)("SQ"), toUTFz!(char *)("SN"), toUTFz!(char *)(key));
    }
    
    // add new SQ records with lengths from previous records
    auto lengths = bamr.header.targetLengths;
    foreach (i, contig; bamr.header.targetNames)
    {
        if(!(contig in mapping)) continue;
        newHeader.addLine(RecordType.SQ, "SN", mapping[contig], "LN", lengths[i].to!string);
    }
    // add PG record
    newHeader.addLine(RecordType.PG, "ID", "recontig", "PN", "recontig", "VN", "1.0.0", "CL", argStr);

    // set writer
    static if(outputBam) auto bamw =  SAMWriter("-",newHeader,SAMWriterTypes.BAM);
    else auto bamw =  SAMWriter("-",newHeader,SAMWriterTypes.SAM);

    // loop over records and convert tid and matetid
    // then write
    auto contigs = bamr.header.targetNames;
    foreach (SAMRecord rec; bamr.allRecords)
    {
        // if unmapped check mate and write
        if(!rec.isMapped){
            if(rec.isMateMapped){
                auto newMateContig = mapping[contigs[rec.mateTID]];
                rec.mateTID = bamw.header.targetId(newMateContig);
            }
            bamw.write(rec);

        // if newcontig doesnt exist write to eject sam
        }else if(!(contigs[rec.tid] in mapping)){ 
            ejectedbamw.write(rec);

        // else convert tid and check mate and write
        }else{
            if(rec.isMateMapped){
                auto newMateContig = mapping[contigs[rec.mateTID]];
                rec.mateTID = bamw.header.targetId(newMateContig);
            }
            auto newContig = mapping[contigs[rec.tid]];
            rec.tid = bamw.header.targetId(newContig);
            bamw.write(rec);
        }
    }
}