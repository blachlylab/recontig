module recontig.bed;

import dhtslib.bed;
import recontig : HEADER_MOD, recontigLine;

void recontigBed(string fn, string ejectedfn, string[string] mapping, string fileOut, string argStr)
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto bedr = BedReader(fn);
    string header;
    if(bedr.header == "") header = HEADER_MOD ~ argStr;
    else header = bedr.header ~ "\n" ~ HEADER_MOD ~ argStr;
    auto bedw = BedWriter(fileOut, header);
    auto ejectedbedw = BedWriter(ejectedfn, bedr.header);

    foreach (BedRecord rec; bedr)
    {
        if(!(rec.contig in mapping)) ejectedbedw.write(rec);
        else {
            rec.contig = mapping[rec.contig];
            bedw.write(rec);
        }
    }

}

string recontigBedRecord(string line, string[string] mapping)
{
    return recontigLine(line, 0, mapping);
}