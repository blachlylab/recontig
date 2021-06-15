import std.stdio;
import std.getopt;
import std.conv : to;
import std.array : join;
import std.format : format;
import std.traits : EnumMembers;
import std.algorithm : map;

import dhtslib.bgzf;
import htslib.hts_log;
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
string fileOut = "-";
bool verbose;
bool verbose2;
bool quiet;
int col = 0;
string delimiter;
string HELP =  
"recontig: remap contig names for different bioinformatics file types.
                
usage: recontig [-e ejected.txt] [-o output] [-m mapping.txt | -b build -c conversion] [-f filetype | --col 1 --delimiter ','] <in.file>
Input can be any of the following formats: vcf, bcf, bam, sam, bed, gff
Input can also be a delimited record based file 
Input can be compressed with gzip or bgzf and can be accessed remotely via https or s3 (see htslib for details).
use 'recontig build-help' to check availiable builds                                 
use 'recontig -b build' conversion-help to check availiable conversions for a build
";

int main(string[] args)
{
	auto clstr = args.join(" ");
	auto res = getopt(args, 
			"file-type|f", "Type of file to convert (vcf, bcf, bam, sam, bed, gff)", &type,
			"col", "if converting a generic file you can specify a column", &col,
			"delimiter", "if converting a generic file you can specify a delimiter (default: '\\t')", &delimiter,
			"build|b","Genome build i.e GRCh37 for using dpryan79's files", &build, 
			"conversion|c", "Conversion string i.e UCSC2ensembl for using dpryan79's files", &conversion,
			"mapping|m", "If want to use your own remapping file instead of dpryan79's", &mappingfn,
			"ejected-output|e", "File to write ejected records to (records with unmapped contigs)", &ejectedfn,
			"verbose|v", "print extra information", &verbose,
			"debug", "print extra debug information", &verbose2,
			"quiet|q", "silence warnings", &quiet,
			"output|o", "name of file out (default is - for stdout)", &fileOut,
		);
	hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
	if(quiet) hts_set_log_level(htsLogLevel.HTS_LOG_ERROR);
	if(verbose) hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
	if(verbose2) hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
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
			hts_log_error("recontig","Error: -b was not set. Need valid build to provide availiable conversions");
			return 1;
		}
		bool buildFound;
		ulong buildIdx = 0;
		foreach(i, b;BUILDS){
			if(b == build) buildFound = true, buildIdx = i;
		}
		if(!buildFound){
			hts_log_error("recontig","Error: Please use a valid build: " ~ BUILDS.to!string);
			return 1;
		}
		stderr.writeln("Valid conversions for " ~ build ~ ": " ~ CONVERSIONS[buildIdx].to!string);
		return 0;
	}

	if(args[1] == "make-mapping"){
		if(args.length != 4){
			hts_log_error("recontig","Error: Need two fastas to make mapping file");
			return 1;
		}
		makeMapping(args[2], args[3]);
		return 0;
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
				ejectedfn = "ejected.txt";
		}
	}

	/// load mapping
	/// if a mapping file is provided load it 
	/// else get a dpryan79 file
	string[string] mapping;
	if(mappingfn != "") mapping = getContigMapping(mappingfn);
	else{
		if(build =="" || conversion == ""){
			hts_log_error("recontig","Error: if not using a mapping file you must provide a valid build and conversion.");
			return 1;
		}
		// validate build
		bool buildFound;
		ulong buildIdx = 0;
		foreach(i, b;BUILDS){
			if(b == build) buildFound = true, buildIdx = i;
		}
		if(!buildFound){
			hts_log_error("recontig","Error: Please use a valid build: " ~ BUILDS.to!string);
			return 1;
		}

		//validate conversion
		bool convFound;
		ulong convIdx = 0;
		foreach(i, c;CONVERSIONS[buildIdx]){
			if(c == conversion) convFound = true, convIdx = i;
		}
		if(!convFound){
			hts_log_error("recontig","Error: Please use a valid conversion: " ~ CONVERSIONS[buildIdx].to!string);
			return 1;
		}
		mapping = getDpryan79ContigMapping(build, conversion);
	}
	

	/// lets recontig
	switch(type)
	{
		case InputFileType.vcf:
		case InputFileType.bcf:
			recontigVcf(args[1], ejectedfn, mapping, fileOut, clstr);
			break;
		case InputFileType.gff:
			recontigGff(args[1], ejectedfn, mapping, fileOut, clstr);
			break;
		case InputFileType.bed:
			recontigBed(args[1], ejectedfn, mapping, fileOut, clstr);
			break;
		case InputFileType.bam:
			recontigBam(args[1], ejectedfn, mapping, fileOut, clstr);
			break;
		case InputFileType.sam:
			recontigSam(args[1], ejectedfn, mapping, fileOut, clstr);
			break;
		default:
			if(col--)
				recontigGeneric(args[1], ejectedfn, col, mapping, fileOut, delimiter);
			else{
				hts_log_error("recontig","Error: Filetype must be specified or -c must be used for a generic file type.");
				return 1;
			}
	}
	return 0;
}