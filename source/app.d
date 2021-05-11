import std.stdio;
import std.getopt;
import std.conv : to;

import mapping;
import vcf;

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

	switch(type)
	{
		case InputFileType.vcf:
		case InputFileType.bcf:
			recontigVCF(args, build, conversion);
			break;
		case InputFileType.gff:
		case InputFileType.bed:
		case InputFileType.bam:
		case InputFileType.sam:
			stderr.writeln("Error: Filetype not supported yet.");
			return;
		default:
			stderr.writeln("Error: Filetype not supported.");
			return;
	}

	
}
