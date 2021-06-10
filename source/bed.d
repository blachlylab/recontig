module bed;

import dhtslib.bed;

void recontigBed(string fn, string ejectedfn, string[string] mapping, string argStr = "")
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto bedr = BedReader(fn);
    auto bedw = BedWriter("-", bedr.header);
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