module recontig.vcf;

import std.utf : toUTFz;
import std.format : format;

import dhtslib.vcf;
import htslib.vcf;


void recontigVcf(string fn, string ejectedfn, string[string] mapping, string argStr)
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto vcfr = VCFReader(fn);
	auto ejectedvcfw = VCFWriter(ejectedfn, vcfr.getHeader);

	// get headers
	auto oldHeader = vcfr.getHeader;
	auto newHeader = VCFHeader(bcf_hdr_dup(oldHeader.hdr));
	// clean new header
	bcf_hdr_remove(newHeader.hdr, BCF_HL_CTG, null);
	foreach (ctg; oldHeader.sequences)
	{
		// if no mapping, skip
		if(!(ctg in mapping)) continue;

		// get header record for contig
		auto hdrRec = bcf_hrec_dup(bcf_hdr_get_hrec(oldHeader.hdr, BCF_HL_CTG, toUTFz!(const(char) *)("ID"), toUTFz!(const(char) *)(ctg), null));

		// get ID key and set value to new contig mapping
		auto key = bcf_hrec_find_key(hdrRec, toUTFz!(const(char) *)("ID"));
		if(key == -1) throw new Exception("hdr_hrec_find_key failed");
		auto err = bcf_hrec_set_val(hdrRec, key, toUTFz!(const(char) *)(mapping[ctg]), mapping[ctg].length, 0);
		if(err == -1) throw new Exception("hdr_hrec_set_value failed");

		// add new record to new header
		bcf_hdr_add_hrec(newHeader.hdr, hdrRec);
	}

	// sync header
	bcf_hdr_sync(newHeader.hdr);

	// make vcfwriter and write header
	auto vcfw = VCFWriter("-", &newHeader);
    if(oldHeader.sequences.length == 0)
    {
        foreach (ctg; mapping.byValue)
        {
            if(ctg == "") continue;
            vcfw.addHeaderLineRaw("##contig=<ID=%s>".format(ctg));
        }
    }
	vcfw.addHeaderLineRaw("##source=recontig");
	vcfw.addHeaderLineRaw("##recontigCMD=<ID=recontig,VERSION=1.0.0,CMDLINE=\"%s\">".format(argStr));
	vcfw.writeHeader;
	ejectedvcfw.writeHeader;

	// loop over records from reader
	foreach (rec; vcfr)
	{
		// if chrom not in mapping, skip
		if(!(rec.chrom in mapping)) ejectedvcfw.writeRecord(rec);
		else{
			// get old chrom, set new header, remap, and write
			auto oldchrom = rec.chrom;
			rec.vcfheader = vcfw.getHeader;
			rec.chrom = mapping[oldchrom];
			vcfw.writeRecord(rec);
		}
	}
}