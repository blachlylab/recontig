import std.stdio;
import std.utf : toUTFz;
import std.getopt;
import std.conv : to;
import dhtslib.vcf;
import htslib.vcf;
import mapping;

string build;
string conversion;
string type;

// dhtslib's version of this function isn't correctly bound
extern (C) int bcf_hrec_set_val(bcf_hrec_t *hrec, int i, const(char) *str, size_t len, int is_quoted);

void main(string[] args)
{
	auto res = getopt(args, 
			config.required, "build|b","Genome build i.e GRCh37", &build, 
			config.required, "conversion|c", "Conversion string i.e UCSC2ensembl", &conversion,
			config.required, "file-type|f", "Type of file to convert i.e vcf", &type
		);
	if (res.helpWanted | (args.length < 2))
	{
		defaultGetoptPrinter("recontig: remap contig names for different bioinformatic file types",
				res.options);
		stderr.writeln();
		return;
	}
	
	// validate build
	bool buildFound;
	ulong buildIdx = 0;
	foreach(i, b;BUILDS){
		if(b == build) buildFound = true, buildIdx = i;
	}
	if(!buildFound) throw new Exception("Please use a valid build: " ~ BUILDS.to!string);

	//validate conversion
	bool convFound;
	ulong convIdx = 0;
	foreach(i, c;CONVERSIONS[buildIdx]){
		if(c == conversion) convFound = true, convIdx = i;
	}
	if(!convFound) throw new Exception("Please use a valid conversion: " ~ CONVERSIONS[buildIdx].to!string);

	// get mapping
	auto mapping = getContigMapping(build, conversion);

	// open reader
	auto vcfr = VCFReader(args[1]);

	// get headers
	auto oldHeader = vcfr.getHeader;
	auto newHeader = VCFHeader(bcf_hdr_dup(oldHeader.hdr));
	// clean new header
	bcf_hdr_remove(newHeader.hdr, BCF_HL_CTG, null);
	foreach (ctg; oldHeader.sequences)
	{
		// if no mapping, skip
		if(mapping[ctg] == "") continue;

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
	vcfw.writeHeader;

	// loop over records from reader
	foreach (rec; vcfr)
	{
		// if chrom not in mapping, skip
		if(mapping[rec.chrom] == "") continue;

		// get old chrom, set new header, remap, and write
		auto oldchrom = rec.chrom;
		rec.vcfheader = vcfw.getHeader;
		rec.chrom = mapping[oldchrom];
		vcfw.writeRecord(rec);
	}
}
