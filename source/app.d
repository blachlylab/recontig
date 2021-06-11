import std.stdio;
import std.getopt;
import std.conv : to;
import std.array : join;
import std.format : format;
import std.traits : EnumMembers;
import std.algorithm : map;

import dhtslib.bgzf;
import recontig;

enum InputFileType
{
	None,
	vcf,
	bcf,
	bam,
	sam,
	bed,
	gff
}

string build;
string conversion;
InputFileType type = InputFileType.None;
string ejectedfn;
string mappingfn;
string HELP =  
"recontig: remap contig names for different bioinformatics file types.
                
usage: recontig [-e ejected.txt] [-m mapping.txt | -b build -c conversion] <in.file>

Input can be any of the following formats: vcf, bcf, bam, sam, bed, gff              
Input can be compressed with gzip or bgzf and can be accessed remotely via https or s3 (see htslib for details).

use 'recontig build-help' to check availiable builds                                 
use 'recontig -b build' conversion-help to check availiable conversions for a build
";

int main(string[] args)
{
	auto clstr = args.join(" ");
	auto res = getopt(args, 
			"file-type|f", "Type of file to convert (vcf, bcf, bam, sam, bed, gff)", &type,
			"build|b","Genome build i.e GRCh37 for using dpryan79's files", &build, 
			"conversion|c", "Conversion string i.e UCSC2ensembl for using dpryan79's files", &conversion,
			"mapping|m", "If want to use your own remapping file instead of dpryan79's", &mappingfn,
			"ejected-output|e", "File to write ejected records to (records with unmapped contigs)", &ejectedfn
		);
	
	if (res.helpWanted || (args.length < 2))
	{
		defaultGetoptPrinter(HELP,
				res.options);
		stderr.writeln();
		return 0;
	}

	if(args[1] == "build-help"){
		stderr.writeln("Valid builds: " ~ BUILDS.to!string);
		return 0;
	}
	
	if(args[1] == "conversion-help"){
		if(build == ""){
			stderr.writeln("Error: -b was not set. Need valid build to provide availiable conversions");
			return 1;
		}
		bool buildFound;
		ulong buildIdx = 0;
		foreach(i, b;BUILDS){
			if(b == build) buildFound = true, buildIdx = i;
		}
		if(!buildFound){
			stderr.writeln("Error: Please use a valid build: " ~ BUILDS.to!string);
			return 1;
		}
		stderr.writeln("Valid conversions for " ~ build ~ ": " ~ CONVERSIONS[buildIdx].to!string);
		return 0;
	}

	if(args[1] == "make-mapping"){
		if(args.length != 4){
			stderr.writeln("Error: Need two fastas to make mapping file");
			return 1;
		}
		makeMapping(args[2], args[3]);
		return 0;
	}

	if (type == InputFileType.None)
	{
		defaultGetoptPrinter("Error: -f or --file-type is required.\n" ~ HELP,
				res.options);
		return 1;
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
				return 1;
		}
	}

	/// if a mapping file is provided load it 
	/// else get a dpryan79 file
	BGZFile mappingFile;
	if(mappingfn != "") mappingFile = BGZFile(mappingfn);
	else{
		if(build =="" || conversion == ""){
			stderr.writeln("Error: if not using a mapping file you must provide a valid build and conversion.");
			return 1;
		}
		// validate build
		bool buildFound;
		ulong buildIdx = 0;
		foreach(i, b;BUILDS){
			if(b == build) buildFound = true, buildIdx = i;
		}
		if(!buildFound){
			stderr.writeln("Error: Please use a valid build: " ~ BUILDS.to!string);
			return 1;
		}

		//validate conversion
		bool convFound;
		ulong convIdx = 0;
		foreach(i, c;CONVERSIONS[buildIdx]){
			if(c == conversion) convFound = true, convIdx = i;
		}
		if(!convFound){
			stderr.writeln("Error: Please use a valid conversion: " ~ CONVERSIONS[buildIdx].to!string);
			return 1;
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
			break;
		default:
			stderr.writeln("Error: Filetype not supported.");
			return 1;
	}
	return 0;
}
