module tests;

import std.process;
import std.path;

string runRecontig(string[] args)
{
    auto binary = "../recontig";
    auto ps = execute([binary] ~ args, null, Config(Config.Flags.stderrPassThrough));
    assert(ps.status == 0);
    return ps.output;
}

/// test binary
unittest
{
	runRecontig([]);
}