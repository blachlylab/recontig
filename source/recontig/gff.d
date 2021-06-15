module recontig.gff;

import dhtslib.gff;
import recontig : HEADER_MOD, recontigLine;

void recontigGff(string fn, string ejectedfn, string[string] mapping, string argStr)
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto gffr = GFF2Reader(fn);
    string header;
    if(gffr.header == "") header = HEADER_MOD ~ argStr;
    else header = gffr.header ~ "\n" ~ HEADER_MOD ~ argStr;
    auto gffw = GFF2Writer("-", header);
    auto ejectedgffw = GFF2Writer(ejectedfn, gffr.header);

    foreach (GFF2Record rec; gffr)
    {
        if(!(rec.contig in mapping)) ejectedgffw.write(rec);
        else {
            rec.contig = mapping[rec.contig];
            gffw.write(rec);
        }
    }

}

string recontigGFFRecord(string line, string[string] mapping)
{
    return recontigLine(line, 0, mapping);
}