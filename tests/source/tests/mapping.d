module tests.mapping;

import std.stdio;
import std.algorithm: map, each;
import std.array: split;
import tests;

/// test binary
unittest
{
	auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test1.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    writeln(res);
    assert(mapping["chrA"]     == "A_A");
    assert(mapping["chrB"]     == "B_A");
    assert(mapping["chrC"]     == "C_A");
    assert(mapping["chrD"]     == "D_A");
    assert(mapping["chrBC"]   == "BC_A");
    assert(mapping["chrBD"]   == "BD_A");
    assert(mapping["chrBCD"] == "BCD_A");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test2.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_B");
    assert(mapping["chrB"]     == "B_B");
    assert(mapping["chrC"]     == "C_B");
    assert(mapping["chrD"]     == "D_B");
    // assert(mapping["chrBC"]   == "BC_B");
    assert(mapping["chrBD"]   == "BD_B");
    assert(mapping["chrBCD"] == "BCD_B");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test3.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_C");
    assert(mapping["chrB"]     == "B_C");
    assert(mapping["chrC"]     == "C_C");
    assert(mapping["chrD"]     == "D_C");
    // assert(mapping["chrBC"]   == "BC_C");
    assert(mapping["chrBD"]   == "BD_C");
    assert(mapping["chrBCD"] == "BCD_C");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test4.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_D");
    assert(mapping["chrB"]     == "B_D");
    assert(mapping["chrC"]     == "C_D");
    assert(mapping["chrD"]     == "D_D");
    assert(mapping["chrBC"]   == "BC_D");
    assert(mapping["chrBD"]   == "BD_D");
    assert(mapping["chrBCD"] == "BCD_D");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test5.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    assert(mapping["chrA"]     == "A_BC");
    assert(mapping["chrB"]     == "B_BC");
    assert(mapping["chrC"]     == "C_BC");
    assert(mapping["chrD"]     == "D_BC");
    assert(mapping["chrBC"]   == "BC_BC");
    assert(mapping["chrBD"]   == "BD_BC");
    assert(mapping["chrBCD"] == "BCD_BC");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test6.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_BD");
    assert(mapping["chrB"]     == "B_BD");
    assert(mapping["chrC"]     == "C_BD");
    assert(mapping["chrD"]     == "D_BD");
    assert(mapping["chrBC"]   == "BC_BD");
    assert(mapping["chrBD"]   == "BD_BD");
    assert(mapping["chrBCD"] == "BCD_BD");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test7.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_CD");
    assert(mapping["chrB"]     == "B_CD");
    assert(mapping["chrC"]     == "C_CD");
    assert(mapping["chrD"]     == "D_CD");
    assert(mapping["chrBC"]   == "BC_CD");
    assert(mapping["chrBD"]   == "BD_CD");
    assert(mapping["chrBCD"] == "BCD_CD");
}

unittest
{
    auto res = runRecontig(["make-mapping", "--enforce-md5sums","data/baseline.fa","data/test8.fa"]);
    string[string] mapping;
    writeln(res);
    if(res.split("\n").length > 0){
        res.split("\n")[0..$-1]
            .map!(x => x.split("\t"))
            .each!(x=>mapping[x[0]] = x[1]);
    
    }
    
    assert(mapping["chrA"]     == "A_BCD");
    assert(mapping["chrB"]     == "B_BCD");
    assert(mapping["chrC"]     == "C_BCD");
    assert(mapping["chrD"]     == "D_BCD");
    assert(mapping["chrBC"]   == "BC_BCD");
    assert(mapping["chrBD"]   == "BD_BCD");
    assert(mapping["chrBCD"] == "BCD_BCD");
}