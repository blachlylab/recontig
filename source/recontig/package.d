module recontig;

public import recontig.vcf;
public import recontig.bed;
public import recontig.gff;
public import recontig.bam;
public import recontig.mapping;

import std.array : split, join;
import std.format : format;

import htslib.hts_log;

string recontigLine(string line, int contigCol, string[string] mapping, string delimiter="\t")
{
    auto fields = line.split(delimiter);
    auto contig = fields[contigCol];
    if(contig in mapping){
        fields[contigCol] = mapping[contig];
    }else{
        hts_log_warning("recontig","Contig %s not in mapping");
        return "";
    }
    return fields.join(delimiter);
}

import pyd.pyd;

enum HEADER_MOD = "# contig names remapped with recontig version 1.0.0. cmd: ";

extern(C) void PydMain() {
    def!(recontigVcf)();
    def!(recontigVcfHeader)();
    def!(recontigVcfRecord)();
    def!(recontigBed)();
    def!(recontigBedRecord)();
    def!(recontigGff)();
    def!(recontigGffRecord)();
    def!(recontigBam)();
    def!(recontigSam)();
    def!(recontigSamHeader)();
    def!(recontigSamRecord)();
    def!(getContigMapping)();
    def!(getDpryan79ContigMapping)();
	module_init();
}
