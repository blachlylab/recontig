module recontig;

public import recontig.vcf;
public import recontig.bed;
public import recontig.gff;
public import recontig.bam;
public import recontig.mapping;
public import recontig._version;

import std.stdio;
import std.array : split, join;
import std.format : format;
import std.algorithm : startsWith;

import htslib.hts_log;
import dhtslib.bgzf;

void recontigGeneric(string fn, string ejectedfn, int contigCol, string[string] mapping, string fileOut, string delimiter="\t", string commentline = "#")
{
    auto f = BGZFile(fn);
    File output;
    if(fileOut == "-"){
        output = stdout;
    }else{
        output = File(fileOut, "w");
    }
    auto ejected = File(ejectedfn, "w");
    foreach(line; f.byLineCopy){
        if(line.startsWith(commentline)){
            output.writeln(line);
            continue;
        }
        auto fields = line.split(delimiter);
        auto contig = fields[contigCol];
        if(contig in mapping){
            fields[contigCol] = mapping[contig];
            output.writeln(fields.join(delimiter));
        }else{
            ejected.writeln(line);
        }
    }
}

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

enum HEADER_MOD = "# contig names remapped with recontig version "~VERSION~". cmd: ";
version(usepyd){
    import pyd.pyd;
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

}