module recontig.mapping.seq;

import std.algorithm : map, sort;
import std.array : array;
import std.traits : isSomeString;
import std.range : ElementType;
import core.bitop : popcnt;

import dhtslib.coordinates;
import htslib.hts;

/// used along with seq_nt16_int & seq_nt16_table to convert
/// all lowercase nucleotides and non ACGTN nucleotides
/// to N
char[5] iupacToACGTN = ['A','C','G','T','N'];

/// seq_nt16_str with support for soft-masked nucleotides
const(char)[32] seq_nt16_str = ['=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N','=','a','c','m','g','r','s','v','t','w','y','h','k','d','b','n'];

/// seq_nt16_int with support for soft-masked nucleotides
const(int)[32] seq_nt16_int = [ 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 , 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 ];

/// seq_nt16_table with support for soft-masked nucleotides
const(ubyte)[256] seq_nt16_table = [
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    1,     2,     4,     8,     15,    15,    15,    15,    15,    15,    15,    15,    15,     0,    15,    15,
    15,     1,    14,     2,    13,    15,    15,     4,    11,    15,    15,    12,    15,     3,    15,    15,
    15,    15,     5,     6,     8,    15,     7,     9,    15,    10,    15,    15,    15,    15,    15,    15,
    15,  1+16, 14+16,  2+16, 13+16,    15,    15,  4+16, 11+16,    15,    15, 12+16,    15,  3+16, 15+16,    15,
    15,    15,  5+16,  6+16,  8+16,    15,  7+16,  9+16,    15, 10+16,    15,    15,    15,    15,    15,    15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
];

const(ubyte)[256] seq_table_deg_2_N = [
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    1,     2,     4,     8,     15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,     1,    15,     2,    15,    15,    15,     4,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    15,     8,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,    15,
    15,  1+16, 15+16,  2+16, 15+16,    15,    15,  4+16, 15+16,    15,    15, 15+16,    15, 15+16,    15,    15,
    15,    15, 15+16, 15+16,  8+16,    15, 15+16, 15+16,    15, 15+16,    15,    15,    15,    15,    15,    15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
];

const(ubyte)[256] seq_table_soft_2_N = [
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15,15,14,15, 13,15,15,15, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6, 15,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
];

unittest
{
    assert(seq_nt16_str[seq_nt16_table['G']] == 'G');
    assert(seq_nt16_str[seq_nt16_table['A']] == 'A');
    assert(seq_nt16_str[seq_nt16_table['T']] == 'T');
    assert(seq_nt16_str[seq_nt16_table['C']] == 'C');

    assert(seq_nt16_str[seq_nt16_table['g']] == 'g');
    assert(seq_nt16_str[seq_nt16_table['a']] == 'a');
    assert(seq_nt16_str[seq_nt16_table['t']] == 't');
    assert(seq_nt16_str[seq_nt16_table['c']] == 'c');

    assert(seq_nt16_str[seq_nt16_table['c']] == 'c');
    assert(seq_nt16_str[seq_nt16_table['g']] == 'g');
    assert(seq_nt16_str[seq_nt16_table['a']] == 'a');
    assert(seq_nt16_str[seq_nt16_table['t']] == 't');
    assert(seq_nt16_str[seq_nt16_table['=']] == '=');
    assert(seq_nt16_str[seq_nt16_table['m']] == 'm');
    assert(seq_nt16_str[seq_nt16_table['r']] == 'r');
    assert(seq_nt16_str[seq_nt16_table['s']] == 's');
    assert(seq_nt16_str[seq_nt16_table['v']] == 'v');
    assert(seq_nt16_str[seq_nt16_table['w']] == 'w');
    assert(seq_nt16_str[seq_nt16_table['y']] == 'y');
    assert(seq_nt16_str[seq_nt16_table['h']] == 'h');
    assert(seq_nt16_str[seq_nt16_table['k']] == 'k');
    assert(seq_nt16_str[seq_nt16_table['d']] == 'd');
    assert(seq_nt16_str[seq_nt16_table['b']] == 'b');
    assert(seq_nt16_str[seq_nt16_table['n']] == 'n');

    assert(seq_nt16_str[seq_nt16_table['=']] == '=');
    assert(seq_nt16_str[seq_nt16_table['M']] == 'M');
    assert(seq_nt16_str[seq_nt16_table['R']] == 'R');
    assert(seq_nt16_str[seq_nt16_table['S']] == 'S');
    assert(seq_nt16_str[seq_nt16_table['V']] == 'V');
    assert(seq_nt16_str[seq_nt16_table['W']] == 'W');
    assert(seq_nt16_str[seq_nt16_table['Y']] == 'Y');
    assert(seq_nt16_str[seq_nt16_table['H']] == 'H');
    assert(seq_nt16_str[seq_nt16_table['K']] == 'K');
    assert(seq_nt16_str[seq_nt16_table['D']] == 'D');
    assert(seq_nt16_str[seq_nt16_table['B']] == 'B');
    assert(seq_nt16_str[seq_nt16_table['N']] == 'N');
}

pragma(inline, true){ 
    bool isHardMasked(T)(T x){
        return !cast(bool)((seq_nt16_table[x] & 15) ^ 15);
    }

    bool isSoftMasked(T)(T x) {
        return cast(bool)(seq_nt16_table[x] & 16);
    }

    bool isDegenerate(T)(T x) {
        return cast(bool)((popcnt(seq_nt16_table[x] & 15) ^ 1) && !isHardMasked(x)) ; 
    }

    /// Convert soft-masked lowercase nucleotides to upper case
    char softMaskToUpper(T)(T c)
    {
        return seq_nt16_str[seq_nt16_table[cast(char)c] & 15];
    }

    /// Convert soft-masked lowercase nucleotides to hard-masked Ns
    char softMaskToHardMask(T)(T c)
    {
        return seq_nt16_str[seq_table_soft_2_N[cast(char)c]];
    }

    /// Convert soft-masked lowercase nucleotides to hard-masked Ns
    char degenerateToHardMask(T)(T c)
    {
        return seq_nt16_str[seq_table_deg_2_N[cast(char)c]];
    }
}

unittest
{
    import std.stdio;
    assert('c'.isSoftMasked);
    assert('g'.isSoftMasked);
    assert('a'.isSoftMasked);
    assert('t'.isSoftMasked);
    assert('m'.isSoftMasked);
    assert('r'.isSoftMasked);
    assert('s'.isSoftMasked);
    assert('v'.isSoftMasked);
    assert('w'.isSoftMasked);
    assert('y'.isSoftMasked);
    assert('h'.isSoftMasked);
    assert('k'.isSoftMasked);
    assert('d'.isSoftMasked);
    assert('b'.isSoftMasked);
    assert('n'.isSoftMasked);
    assert(!'='.isSoftMasked);

    assert('N'.isHardMasked);
    assert('n'.isHardMasked);

    assert(!'c'.isHardMasked);
    assert(!'g'.isHardMasked);
    assert(!'a'.isHardMasked);
    assert(!'t'.isHardMasked);
    assert(!'='.isHardMasked);
    assert(!'m'.isHardMasked);
    assert(!'r'.isHardMasked);
    assert(!'s'.isHardMasked);
    assert(!'v'.isHardMasked);
    assert(!'w'.isHardMasked);
    assert(!'y'.isHardMasked);
    assert(!'h'.isHardMasked);
    assert(!'k'.isHardMasked);
    assert(!'d'.isHardMasked);
    assert(!'b'.isHardMasked);
    assert(!'='.isHardMasked);
    assert(!'M'.isHardMasked);
    assert(!'R'.isHardMasked);
    assert(!'S'.isHardMasked);
    assert(!'V'.isHardMasked);
    assert(!'W'.isHardMasked);
    assert(!'Y'.isHardMasked);
    assert(!'H'.isHardMasked);
    assert(!'K'.isHardMasked);
    assert(!'D'.isHardMasked);
    assert(!'B'.isHardMasked);
    assert(!'C'.isHardMasked);
    assert(!'G'.isHardMasked);
    assert(!'A'.isHardMasked);
    assert(!'T'.isHardMasked);

    assert('='.isDegenerate);
    assert('M'.isDegenerate);
    assert('R'.isDegenerate);
    assert('S'.isDegenerate);
    assert('V'.isDegenerate);
    assert('W'.isDegenerate);
    assert('Y'.isDegenerate);
    assert('H'.isDegenerate);
    assert('K'.isDegenerate);
    assert('D'.isDegenerate);
    assert('B'.isDegenerate);
    assert(!'N'.isDegenerate);

    assert('m'.isDegenerate);
    assert('r'.isDegenerate);
    assert('s'.isDegenerate);
    assert('v'.isDegenerate);
    assert('w'.isDegenerate);
    assert('y'.isDegenerate);
    assert('h'.isDegenerate);
    assert('k'.isDegenerate);
    assert('d'.isDegenerate);
    assert('b'.isDegenerate);
    assert(!'n'.isDegenerate);

    assert(!'c'.isDegenerate);
    assert(!'g'.isDegenerate);
    assert(!'a'.isDegenerate);
    assert(!'t'.isDegenerate);

    assert(!'C'.isDegenerate);
    assert(!'G'.isDegenerate);
    assert(!'A'.isDegenerate);
    assert(!'T'.isDegenerate);

    assert(!'C'.isSoftMasked);
    assert(!'G'.isSoftMasked);
    assert(!'A'.isSoftMasked);
    assert(!'T'.isSoftMasked);    
    
    assert('c'.softMaskToUpper == 'C');
    assert('g'.softMaskToUpper == 'G');
    assert('a'.softMaskToUpper == 'A');
    assert('t'.softMaskToUpper == 'T');

    assert('c'.softMaskToHardMask == 'N');
    assert('g'.softMaskToHardMask == 'N');
    assert('a'.softMaskToHardMask == 'N');
    assert('t'.softMaskToHardMask == 'N');

    assert('='.degenerateToHardMask == 'N');
    assert('M'.degenerateToHardMask == 'N');
    assert('R'.degenerateToHardMask == 'N');
    assert('S'.degenerateToHardMask == 'N');
    assert('V'.degenerateToHardMask == 'N');
    assert('W'.degenerateToHardMask == 'N');
    assert('Y'.degenerateToHardMask == 'N');
    assert('H'.degenerateToHardMask == 'N');
    assert('K'.degenerateToHardMask == 'N');
    assert('D'.degenerateToHardMask == 'N');
    assert('B'.degenerateToHardMask == 'N');

}

pragma(inline, true){
    /// Clean nucleotide sequences by either:
    /// making all lowercase
    /// remove IUPAC/degenerate nucleotides
    T cleanFasta(alias f, T)(T seq)
    if(isSomeString!T)
    {
        return cast(T) seq.map!(f).array;
    }

    auto convertSoftMaskToUpper(T)(T v){
        return cleanFasta!((ElementType!T x) => x.softMaskToUpper)(v);
    }

    auto convertSoftMaskToHardMask(T)(T v){
        return cleanFasta!((ElementType!T x) => x.softMaskToHardMask)(v);
    }

    auto convertDegenerateToHardMask(T)(T v){
        return cleanFasta!((ElementType!T x) => x.degenerateToHardMask)(v);
    }

    auto convertAllToHardMask(T)(T v){
        return cleanFasta!((ElementType!T x) => x.degenerateToHardMask.softMaskToHardMask)(v);
    }
}
unittest
{
    char[] seq1 = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN".dup;
    assert(seq1.convertSoftMaskToUpper      == "NNNNNNGATCGACTGACTGATCTGA=MRSVWYHKDBGATCGGATCNNNN".dup);
    assert(seq1.convertSoftMaskToHardMask   == "NNNNNNGATCGACTGACTNNNNNNN=MRSVWYHKDBGATCGGATCNNNN");
    assert(seq1.convertDegenerateToHardMask == "NNNNNNGATCGACTGACTgatctgaNNNNNNNNNNNGATCGGATCNNNN");
    assert(seq1.convertAllToHardMask        == "NNNNNNGATCGACTGACTNNNNNNNNNNNNNNNNNNGATCGGATCNNNN");
}

/// template function for recording regions
pragma(inline, true) void getRegions(alias f)(string seq, ZB start, ref ZBHO[] regions, bool * lastRegionFinished)
{
    long i;
    /// Were we in a N-masked region prior?
    if(!(*lastRegionFinished))
    {
        auto last = &regions[$-1];

        while((i < seq.length) && (f(seq[i]))) i++;

        if(i == seq.length){
            *lastRegionFinished = false;
            last.end = start + i;
            return;
        }

        *lastRegionFinished = true;
        last.end = start + i;
    }
    for(; i < seq.length; i++)
    {
        if(f(seq[i])){
            ZBHO reg = ZBHO(start + i, start + (++i));

            while((i < seq.length) && (f(seq[i]))) i++;

            if(i == seq.length){
                *lastRegionFinished = false;
                reg.end = start + i;
                regions ~= reg;
                return;
            }

            *lastRegionFinished = true;
            reg.end = start + i;
            regions ~= reg;
        } 
    }
}


alias getHardMaskedRegions = getRegions!isHardMasked;

unittest
{
    import std.stdio;
    string seq = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN"c;

    ZBHO[] regions;
    bool finished = true;

    getHardMaskedRegions(seq, ZB(0), regions, &finished);

    assert(finished == false);
    assert(regions == [ZBHO(0,6), ZBHO(45,49)]);
    
    regions = [];
    finished = true;
    getHardMaskedRegions(seq, ZB(0), regions, &finished);
    getHardMaskedRegions(seq, ZB(seq.length), regions, &finished);

    assert(finished == false);
    assert(regions == [ZBHO(0,6), ZBHO(45,55), ZBHO(94,98)]);

    regions = [];
    finished = true;
    getHardMaskedRegions(seq, ZB(0), regions, &finished);
    getHardMaskedRegions(seq, ZB(seq.length), regions, &finished);
    getHardMaskedRegions(seq[0..13], ZB(seq.length + seq.length), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(0,6), ZBHO(45,55), ZBHO(94,104)]);

}

alias getSoftMaskedRegions = getRegions!isSoftMasked;

unittest
{
    import std.stdio;
    string seq = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN"c;

    ZBHO[] regions;
    bool finished = true;

    getSoftMaskedRegions(seq, ZB(0), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(18,25)]);

}

alias getDegerateRegions = getRegions!isDegenerate;

unittest
{
    import std.stdio;
    string seq = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN"c;

    ZBHO[] regions;
    bool finished = true;

    getDegerateRegions(seq, ZB(0), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(25, 36)]);

}

unittest
{
    string seq1 = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN".convertDegenerateToHardMask.convertSoftMaskToUpper;
    string seq2 = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN".convertAllToHardMask;

    assert(seq1 ==  "NNNNNNGATCGACTGACTGATCTGANNNNNNNNNNNGATCGGATCNNNN");
    assert(seq2 ==  "NNNNNNGATCGACTGACTNNNNNNNNNNNNNNNNNNGATCGGATCNNNN");

    ZBHO[] regions;
    bool finished = true;

    getHardMaskedRegions(seq1, ZB(0), regions, &finished);

    assert(finished == false);
    assert(regions == [ZBHO(0,6), ZBHO(25, 36), ZBHO(45,49)]);
    
    regions = [];
    finished = true;
    getHardMaskedRegions(seq1, ZB(0), regions, &finished);
    getHardMaskedRegions(seq2, ZB(seq1.length), regions, &finished);

    assert(finished == false);
    assert(regions == [ZBHO(0,6), ZBHO(25, 36), ZBHO(45,55), ZBHO(67, 85), ZBHO(94,98)]);

    regions = [];
    finished = true;
    getHardMaskedRegions(seq1, ZB(0), regions, &finished);
    getHardMaskedRegions(seq2, ZB(seq1.length), regions, &finished);
    getHardMaskedRegions(seq1[0..13], ZB(seq1.length + seq2.length), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(0,6), ZBHO(25, 36), ZBHO(45,55), ZBHO(67, 85), ZBHO(94,104)]);

}

pragma(inline, true) auto hardMaskRegions(char[] seq, ref ZBHO[] regions, ZBHO coords)
{
    while(regions.length > 0 && regions[0].isOverlap(coords))
    {
        auto reg = regions[0];

        auto overlap = reg & coords;

        overlap = overlap.offset(-coords.start);
        
        seq[overlap.start .. overlap.end] = 'N';

        if(coords.end < reg.end) break;
        else regions = regions[1..$];
    }
}

unittest
{
    string seq1 = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN".convertDegenerateToHardMask.convertSoftMaskToUpper;
    string seq2 = "NNNNNNGATCGACTGACTgatctga=MRSVWYHKDBGATCGGATCNNNN".convertAllToHardMask;

    assert(seq1 ==  "NNNNNNGATCGACTGACTGATCTGANNNNNNNNNNNGATCGGATCNNNN");
    assert(seq2 ==  "NNNNNNGATCGACTGACTNNNNNNNNNNNNNNNNNNGATCGGATCNNNN");

    ZBHO[] regions;
    bool finished = true;

    regions = [];
    finished = true;
    getHardMaskedRegions(seq2, ZB(0), regions, &finished);

    assert(finished == false);
    assert(regions == [ZBHO(0,6), ZBHO(18, 36), ZBHO(45,49)]);

    auto mod = seq1.dup;
    hardMaskRegions(mod, regions, ZBHO(0, 20));
    assert(mod == "NNNNNNGATCGACTGACTNNTCTGANNNNNNNNNNNGATCGGATCNNNN".dup);
    hardMaskRegions(mod, regions, ZBHO(0, 21));
    assert(mod == "NNNNNNGATCGACTGACTNNNCTGANNNNNNNNNNNGATCGGATCNNNN".dup);
    hardMaskRegions(mod, regions, ZBHO(0, mod.length));
    assert(mod == seq2.dup);
}

pragma(inline, true) auto convertDegenerateRegions(char[] seq, char[] seq2, ref ZBHO[] regions, ZBHO coords)
{
    while(regions.length > 0 && regions[0].isOverlap(coords))
    {
        auto reg = regions[0];

        auto overlap = reg & coords;

        overlap = overlap.offset(-coords.start);

        ubyte[] bytes1 = seq[overlap.start .. overlap.end].map!(x => cast(ubyte) seq_nt16_table[x]).array;
        ubyte[] bytes2 = seq2[overlap.start .. overlap.end].map!(x => cast(ubyte) seq_nt16_table[x]).array;

        bytes1[] = bytes1[] | bytes2[];

        seq[overlap.start .. overlap.end] = bytes1.map!(x => seq_nt16_str[x]).array;

        if(coords.end < reg.end) break;
        else regions = regions[1..$];
    }
}

unittest
{
    string seq1 = "NNNNNNGATCGACTGACTgatctgaWSMKRYBDHVWSMKRGATCGGATCNNNN";
    string seq2 = "NNNNNNGATCGACTGACTgatctgaACATGCCAAAACATGGATCGGATCNNNN";

    ZBHO[] regions;
    bool finished = true;

    regions = [];
    finished = true;
    getDegerateRegions(seq1, ZB(0), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(25, 40)]);

    
    auto mod = seq1.dup;
    convertDegenerateRegions(mod, seq2.dup, regions, ZBHO(0, seq1.length));
    assert(mod == "NNNNNNGATCGACTGACTgatctgaWSMKRYBDHVWSMKRGATCGGATCNNNN".dup);
}

pragma(inline, true) auto convertSoftMaskedRegions(char[] seq, ref ZBHO[] regions, ZBHO coords)
{
    while(regions.length > 0 && regions[0].isOverlap(coords))
    {
        auto reg = regions[0];

        auto overlap = reg & coords;

        overlap = overlap.offset(-coords.start);

        seq[overlap.start .. overlap.end] = seq[overlap.start .. overlap.end].convertSoftMaskToUpper;

        if(coords.end < reg.end) break;
        else regions = regions[1..$];
    }
}

unittest
{
    string seq1 = "NNNNNNGATCGACTGACTgatctgaACATGCCAAAACATGGATCGGATCNNNN";

    ZBHO[] regions;
    bool finished = true;

    regions = [];
    finished = true;
    getSoftMaskedRegions(seq1, ZB(0), regions, &finished);

    assert(finished == true);
    assert(regions == [ZBHO(18, 25)]);

    
    auto mod = seq1.dup;
    convertSoftMaskedRegions(mod, regions, ZBHO(0, seq1.length));
    assert(mod == "NNNNNNGATCGACTGACTGATCTGAACATGCCAAAACATGGATCGGATCNNNN".dup);
}

/// get a union of regions
ZBHO[] unionRegions(ZBHO[] regs1, ZBHO[] regs2)
{
    ZBHO[] ret;
    auto regs = (regs1 ~ regs2).sort!"a.start < b.start";
    if(regs.empty) return [];
    ZBHO last = regs.front;
    regs.popFront;
    foreach(reg; regs)
    {
        if(last.isOverlap(reg)) last = last | reg;
        else{
            ret ~= last;
            last = reg;
        }
    }
    if(ret.length == 0) ret ~= last;
    else if(last != ret[$-1]) ret ~= last;
    return ret;
}

unittest
{

    ZBHO[] a = [ZBHO(0, 2), ZBHO(3, 5), ZBHO(7, 10)];
    ZBHO[] b = [ZBHO(1, 4), ZBHO(6, 8), ZBHO(7, 10), ZBHO(13, 15)];
    import std.stdio;
    assert(unionRegions(a, b) == [ZBHO(0, 5), ZBHO(6, 10), ZBHO(13, 15)]);
}