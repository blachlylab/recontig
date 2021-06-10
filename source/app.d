import std.stdio;
import std.getopt;
import std.conv : to;
import std.array : join;

import dhtslib.bgzf;
import recontig;

enum InputFileType
{
	vcf,
	bam,
	sam,
	bcf,
	bed,
	gff
}

string build;
string conversion;
InputFileType type;
string ejectedfn;
string mappingfn;

void main(string[] args)
{
	auto clstr = args.join(" ");
	auto res = getopt(args, 
			config.required, "file-type|f", "Type of file to convert i.e vcf", &type,
			"build|b","Genome build i.e GRCh37 for using dpryan79's files", &build, 
			"conversion|c", "Conversion string i.e UCSC2ensembl for using dpryan79's files", &conversion,
			"mapping", "If want to use your own remapping file instead of dpryan79's", &mappingfn,
			"ejected-output", "File to write ejected records to (records with unmapped contigs)", &ejectedfn
		);
	if (res.helpWanted | (args.length < 2))
	{
		defaultGetoptPrinter("recontig: remap contig names for different bioinformatics file types",
				res.options);
		stderr.writeln();
		return;
	}
	
	/// if no ejected filename provided we will select a default name
	if(ejectedfn == "")
	{
		switch(type)
		{
			case InputFileType.vcf:
			case InputFileType.bcf:
				ejectedfn = "ejected.vcf";
				break;
			case InputFileType.gff:
				ejectedfn = "ejected.txt";
				break;
			case InputFileType.bed:
				ejectedfn = "ejected.bed";
				break;
			case InputFileType.bam:
			case InputFileType.sam:
				ejectedfn = "ejected.sam";
				break;
			default:
				stderr.writeln("Error: Filetype not supported.");
				return;
		}
	}

	/// if a mapping file is provided load it 
	/// else get a dpryan79 file
	BGZFile mappingFile;
	if(mappingfn != "") mappingFile = BGZFile(mappingfn);
	else{
		if(build =="" || conversion == ""){
			stderr.writeln("Error: if not using a mapping file you must provide a valid build and conversion.");
			return;
		}
		// validate build
		bool buildFound;
		ulong buildIdx = 0;
		foreach(i, b;BUILDS){
			if(b == build) buildFound = true, buildIdx = i;
		}
		if(!buildFound){
			stderr.writeln("Error: Please use a valid build: " ~ BUILDS.to!string);
			return;
		}

		//validate conversion
		bool convFound;
		ulong convIdx = 0;
		foreach(i, c;CONVERSIONS[buildIdx]){
			if(c == conversion) convFound = true, convIdx = i;
		}
		if(!convFound){
			stderr.writeln("Error: Please use a valid conversion: " ~ CONVERSIONS[buildIdx].to!string);
			return;
		}
		mappingFile = getDpryan79ContigMappingFile(build, conversion);
	}

	/// load mapping
	auto mapping = getContigMapping(mappingFile);

	/// lets recontig
	switch(type)
	{
		case InputFileType.vcf:
		case InputFileType.bcf:
			recontigVcf(args[1], ejectedfn, mapping, clstr);
			break;
		case InputFileType.gff:
			recontigGff(args[1], ejectedfn, mapping, clstr);
			break;
		case InputFileType.bed:
			recontigBed(args[1], ejectedfn, mapping, clstr);
			break;
		case InputFileType.bam:
			recontigBam(args[1], ejectedfn, mapping, clstr);
			break;
		case InputFileType.sam:
			recontigSam(args[1], ejectedfn, mapping, clstr);
			return;
		default:
			stderr.writeln("Error: Filetype not supported.");
			return;
	}
}
