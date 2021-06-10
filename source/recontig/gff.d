module gff;

import dhtslib.gff;

void recontigGFF(string fn, string ejectedfn, string[string] mapping, string argStr = "")
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto gffr = GFF2Reader(fn);
    auto gffw = GFF2Writer("-", gffr.header);
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