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

/// build string
string build;

/// conversion string
string conversion;

/// file type
InputFileType type = InputFileType.None;

/// ejected file name
string ejectedfn;

/// mapping file name
string mappingfn;

/// output file name
string fileOut = "-";

/// verbose level 1
bool verbose;

/// verbose level 2
bool verbose2;

/// verbose level 0
bool quiet;

/// column to recontig
int col = 0;

/// delimiter for generic files
string delimiter;

/// comment start for generic files
string comment = "#";

/// help string
string SUBHELP =  
"recontig: convert contig names for different bioinformatics file types.

Subcommands:
build-help          check availiable builds from dpryan79's github
conversion-help     check availiable conversions for a specified build from dpryan79's github
convert             convert a file from one naming convention to another
make-mapping        make a contig conversion file from two fasta files
";

/// help string
string convHELP =  
"recontig conversion-help: check availiable conversions for a specified build from dpryan79's github

usage: recontig conversion-help -b build 
";

/// help string
string MAPHELP =  
"recontig make-mapping: make a contig conversion file from two fasta files
Makes a mapping file for two fasta files that have been faidx'd (samtools faidx)
Fastas can be compressed with bgzf and can be accessed remotely via https or s3 (see htslib for details).

usage: recontig make-mapping [-o output] <from.fa> <to.fa>
";

/// help string
string CONVERTHELP =  
"recontig convert: remap contig names for different bioinformatics file types.

usage: recontig convert [-e ejected.txt] [-o output] [-m mapping.txt | -b build -c conversion] [-f filetype | --col 1 --delimiter ','] <in.file>

Input can be any of the following formats: vcf, bcf, bam, sam, bed, gff
Input can also be a delimited record based file 
Input can be compressed with gzip or bgzf and can be accessed remotely via https or s3 (see htslib for details).
";

int main(string[] args)
{

	if(args.length <= 1){

		stderr.writeln(SUBHELP);
		return 0;

	}else if(args[1] == "build-help"){

		stderr.writeln("Valid builds: " ~ BUILDS.to!string);
		return 0;

	}else if(args[1] == "conversion-help"){

		return conversionHelp(args);

	}else if(args[1] == "make-mapping"){

		return makeMappingRun(args);
	}else if(args[1] == "convert"){

		return convert(args);
	}else{
		stderr.writefln("Invalid subcommand: %s", args[1]);
		stderr.writeln();
		stderr.writeln(SUBHELP);
		return 1;
	}

}

int conversionHelp(string[] args)
{
	args = args[1..$];
	auto res = getopt(args, 
			"build|b","Genome build i.e GRCh37 for using dpryan79's files", &build, 
			"quiet|q", "silence warnings", &quiet,
			"verbose|v", "print extra information", &verbose,
			"debug", "print extra debug information", &verbose2,
		);
	hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
	if(quiet) hts_set_log_level(htsLogLevel.HTS_LOG_ERROR);
	if(verbose) hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
	if(verbose2) hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
	if (res.helpWanted)
	{
		defaultGetoptPrinter(convHELP,
				res.options);
		stderr.writeln();
		return 0;
	}
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

int makeMappingRun(string[] args)
{
	args = args[1..$];
	auto res = getopt(args, 
			"output|o", "name of file out (default is - for stdout)", &fileOut,
			"quiet|q", "silence warnings", &quiet,
			"verbose|v", "print extra information", &verbose,
			"debug", "print extra debug information", &verbose2,
		);
	hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
	if(quiet) hts_set_log_level(htsLogLevel.HTS_LOG_ERROR);
	if(verbose) hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
	if(verbose2) hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
	if (res.helpWanted || (args.length < 2))
	{
		defaultGetoptPrinter(MAPHELP,
				res.options);
		stderr.writeln();
		return 0;
	}
	if(args.length != 3){
		hts_log_error("recontig","Error: Need two fastas to make mapping file");
		return 1;
	}
	makeMapping(args[1], args[2], fileOut);
	return 0;
}

int convert(string[] args){
	auto clstr = args.join(" ");
	args = args[1..$];
	auto res = getopt(args, 
			"build|b","Genome build i.e GRCh37 for using dpryan79's files", &build, 
			"conversion|c", "Conversion string i.e UCSC2ensembl for using dpryan79's files", &conversion,
			"ejected-output|e", "File to write ejected records to (records with unmapped contigs)", &ejectedfn,
			"file-type|f", "Type of file to convert (vcf, bcf, bam, sam, bed, gff)", &type,
			"mapping|m", "If want to use your own remapping file instead of dpryan79's", &mappingfn,
			"output|o", "name of file out (default is - for stdout)", &fileOut,
			"quiet|q", "silence warnings", &quiet,
			"verbose|v", "print extra information", &verbose,
			"col", "if converting a generic file you can specify a column", &col,
			"comment", "if converting a generic file you can specify what a comment line starts with (default: '#')", &comment,
			"debug", "print extra debug information", &verbose2,
			"delimiter", "if converting a generic file you can specify a delimiter (default: '\\t')", &delimiter,
		);
	hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
	if(quiet) hts_set_log_level(htsLogLevel.HTS_LOG_ERROR);
	if(verbose) hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
	if(verbose2) hts_set_log_level(htsLogLevel.HTS_LOG_DEBUG);
	if (res.helpWanted || (args.length < 2))
	{
		defaultGetoptPrinter(CONVERTHELP,
				res.options);
		stderr.writeln();
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
				recontigGeneric(args[1], ejectedfn, col, mapping, fileOut, delimiter, comment);
			else{
				hts_log_error("recontig","Error: Filetype must be specified or -c must be used for a generic file type.");
				return 1;
			}
	}
	return 0;
} 