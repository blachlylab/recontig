module recontig.vcf;

import std.utf : toUTFz;
import std.string : toStringz, fromStringz;
import std.format : format;

import dhtslib.vcf;
import htslib.vcf;
import htslib.hts_log;
import htslib.kstring;


void recontigVcf(string fn, string ejectedfn, string[string] mapping, string fileOut, string argStr)
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
	bool addedOne;
	foreach (ctg; oldHeader.sequences)
	{
		// if no mapping, skip
		if(!(ctg in mapping)) {
			hts_log_warning("recontig","contig %s not present in mapping".format(ctg));
			continue;
		}

		// get header record for contig
		auto hdrRec = bcf_hrec_dup(bcf_hdr_get_hrec(oldHeader.hdr, BCF_HL_CTG, toUTFz!(const(char) *)("ID"), toUTFz!(const(char) *)(ctg), null));

		// get ID key and set value to new contig mapping
		auto key = bcf_hrec_find_key(hdrRec, toUTFz!(const(char) *)("ID"));
		if(key == -1) throw new Exception("hdr_hrec_find_key failed");
		auto err = bcf_hrec_set_val(hdrRec, key, toUTFz!(const(char) *)(mapping[ctg]), mapping[ctg].length, 0);
		if(err == -1) throw new Exception("hdr_hrec_set_value failed");

		// add new record to new header
		bcf_hdr_add_hrec(newHeader.hdr, hdrRec);
		addedOne = true;
	}

	// if there are contig lines and we remapped none
	// error
	if(!addedOne && oldHeader.sequences.length != 0){
		hts_log_error("recontig","No existing contigs are able to be mapped (Are you using the correct mapping file?)");
		return;
	}
	// sync header
	bcf_hdr_sync(newHeader.hdr);

	// make vcfwriter and write header
	auto vcfw = VCFWriter(fileOut, newHeader);
	// if old header had no seq files
    if(oldHeader.sequences.length == 0)
    {
		hts_log_debug("recontig","no existing contig lines in header, adding");
		hts_log_warning("recontig","your vcf has no existing contig lines");
		hts_log_warning("recontig","recontig cannot check if your mapping file is valid for this vcf");
        foreach (ctg; mapping.byKeyValue)
        {
            vcfw.addHeaderLineRaw("##contig=<ID=%s>".format(ctg.value));
			ejectedvcfw.addHeaderLineRaw("##contig=<ID=%s>".format(ctg.key));
        }
    }

	// add cmd lines
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

string recontigVcfHeader(string fn, string[string] mapping, string argStr)
{
    // get mapping
	// auto mapping = getContigMapping(build, conversion);

	// open reader
	auto vcfr = VCFReader(fn);

	// get headers
	auto oldHeader = vcfr.getHeader;
	auto newHeader = VCFHeader(bcf_hdr_dup(oldHeader.hdr));
	// clean new header
	bcf_hdr_remove(newHeader.hdr, BCF_HL_CTG, null);
	bool addedOne;
	foreach (ctg; oldHeader.sequences)
	{
		// if no mapping, skip
		if(!(ctg in mapping)) {
			hts_log_warning("recontig","contig %s not present in mapping".format(ctg));
			continue;
		}

		// get header record for contig
		auto hdrRec = bcf_hrec_dup(bcf_hdr_get_hrec(oldHeader.hdr, BCF_HL_CTG, toUTFz!(const(char) *)("ID"), toUTFz!(const(char) *)(ctg), null));

		// get ID key and set value to new contig mapping
		auto key = bcf_hrec_find_key(hdrRec, toUTFz!(const(char) *)("ID"));
		if(key == -1) throw new Exception("hdr_hrec_find_key failed");
		auto err = bcf_hrec_set_val(hdrRec, key, toUTFz!(const(char) *)(mapping[ctg]), mapping[ctg].length, 0);
		if(err == -1) throw new Exception("hdr_hrec_set_value failed");

		// add new record to new header
		bcf_hdr_add_hrec(newHeader.hdr, hdrRec);
		addedOne = true;
	}

	// if there are contig lines and we remapped none
	// error
	if(!addedOne && oldHeader.sequences.length != 0){
		hts_log_warning("recontig","No existing contigs are able to be mapped (Are you using the correct mapping file?)");
		return "";
	}
	// sync header
	bcf_hdr_sync(newHeader.hdr);

	// make vcfwriter and write header
	auto vcfw = VCFWriter("-", newHeader);
	// if old header had no seq files
    if(oldHeader.sequences.length == 0)
    {
		hts_log_debug("recontig","no existing contig lines in header, adding");
		hts_log_warning("recontig","your vcf has no existing contig lines");
		hts_log_warning("recontig","recontig cannot check if your mapping file is valid for this vcf");
        foreach (ctg; mapping.byKeyValue)
        {
            auto res = bcf_hdr_append(newHeader.hdr, toStringz("##contig=<ID=%s>".format(ctg.value)));
			if(res) hts_log_warning("recontig","could not add header line for contig %s".format(ctg.value));
			bcf_hdr_sync(newHeader.hdr);
        }
    }
	kstring_t ks;
	auto res = bcf_hdr_format(newHeader.hdr, 0, &ks);
	if(res){
		hts_log_warning("recontig","could not format header");
		return "";
	}
	return fromStringz(ks_c_str(&ks)).idup;
}

string recontigVcfRecord(string vcfRec, string header, string[string] mapping)
{
	auto hdr = bcf_hdr_init(toStringz("w"));
	auto res = bcf_hdr_parse(hdr, toUTFz!(char *)(header));
	if(res){
		hts_log_warning("recontig","could parse header string");
		return "";
	}
	auto vcfHdr = VCFHeader(hdr);
	auto rec = VCFRecord(vcfHdr, vcfRec);
	// if chrom not in mapping, skip
	if(!(rec.chrom in mapping)){
		hts_log_warning("recontig", "contig %s not found in mapping".format(rec.chrom));
		return "";
	}
	else{
		// get old chrom, set new header, remap, and write
		auto oldchrom = rec.chrom;
		rec.vcfheader = vcfHdr;
		rec.chrom = mapping[oldchrom];
	}
	return rec.toString();
}