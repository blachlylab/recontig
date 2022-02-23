module recontig.mapping;

import std.stdio;
import std.array : split;
import std.format : format;
import std.conv : to;

import dhtslib.bgzf;
import htslib.hts_log;

public import recontig.mapping.download;
public import recontig.mapping.generate;

/// load a contig mapping file
auto getContigMapping(string fn)
{
    return BGZFile(fn).convertMappingToHashMap;
}

/// load a contig mapping file into a hashmap
auto convertMappingToHashMap(BGZFile file)
{
    string[string] mapping;
    foreach (line; file.byLineCopy)
    {
        auto fields = line.split("\t");
        if(fields[0] == "" || fields[1] == "")
            hts_log_info("recontig", "line %s in mapping file skipped".format(line));
        mapping[fields[0]] = fields[1];
    }
    return mapping;
}