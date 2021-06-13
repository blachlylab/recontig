module recontig;

public import recontig.vcf;
public import recontig.bed;
public import recontig.gff;
public import recontig.bam;
public import recontig.mapping;

import pyd.pyd;

enum HEADER_MOD = "# contig names remapped with recontig version 1.0.0. cmd: ";

extern(C) void PydMain() {
    def!(recontigVcf)();
    def!(recontigBed)();
    def!(recontigGff)();
    def!(recontigBam)();
    def!(getContigMapping)();
    def!(getDpryan79ContigMapping)();
	module_init();
}
