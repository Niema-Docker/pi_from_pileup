// Given a pile-up file from samtools mpileup, compute the pi diversity statistic (https://doi.org/10.1093/ve/vey041)
// Version 1.0.3
// Compile: g++ -O3 -std=c++11 -o pi_from_pileup pi_from_pileup.cpp

// constants
#define DEFAULT_MIN_DEPTH 10

// includes
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// main function
int main(int argc, char* argv[]) {
    // parse user args
    unsigned long MIN_DEPTH;
    if(argc == 3) {
        MIN_DEPTH = stoul(string(argv[2]));
    } else if(argc != 2 || strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"-help") == 0 || strcmp(argv[1],"--help") == 0) {
        cerr << "USAGE: " << argv[0] << " <input_pileup_file> [min_depth=" << DEFAULT_MIN_DEPTH << ']' << endl; exit(1);
    } else {
        MIN_DEPTH = DEFAULT_MIN_DEPTH;
    }

    // open input file
    ifstream infile(argv[1]);
    if(!infile.good()) {
        cerr << "ERROR: Unable to open input file: " << argv[1] << endl; exit(1);
    }

    // temporary variables for parsing input file
    string tmp;              // temporary string
    string line;             // current line
    string chrom;            // current chomosome label
    unsigned long pos;       // current position
    char ref_nuc;            // current reference nucleotide
    unsigned long num_reads; // current number of reads
    string bases;            // current bases

    // temporary variables for computing D_l at any given locus
    unsigned long n_match = 0;      // number of matches
    unsigned long n_mismatch_A = 0; // number of A mismatches
    unsigned long n_mismatch_C = 0; // number of C mismatches
    unsigned long n_mismatch_G = 0; // number of G mismatches
    unsigned long n_mismatch_T = 0; // number of T mismatches
    unsigned long N;                // temporary variable for storing N (total num bases at this locus)
    long double N_Nminus1;          // temporary variable for storing N(N-1)
    long double D_l;                // temporary variable for storing D_l
    long double pi = 0;             // overall output pi (keep adding to it, and then divide at the end)

    // compute pi statistic
    unsigned long L = 0; // total number of positions
    while(getline(infile,line)) {
        // clean up to prep for next line
        n_match = 0;
        n_mismatch_A = 0;
        n_mismatch_C = 0;
        n_mismatch_G = 0;
        n_mismatch_T = 0;

        // parse next line
        istringstream ss(line);
        getline(ss, chrom, '\t');
        getline(ss, tmp, '\t'); pos = stoul(tmp);
        getline(ss, tmp, '\t'); ref_nuc = tmp[0];
        getline(ss, tmp, '\t'); num_reads = stoul(tmp);
        getline(ss, bases, '\t');

        // count nucleotides at this position
        for(unsigned int i = 0; i < bases.size(); ++i) {
            switch(bases[i]) {
                case '.':
                case ',': ++n_match; break;
                case 'A':
                case 'a': ++n_mismatch_A; break;
                case 'C':
                case 'c': ++n_mismatch_C; break;
                case 'G':
                case 'g': ++n_mismatch_G; break;
                case 'T':
                case 't': ++n_mismatch_T; break;
            }
        }
        N = n_match + n_mismatch_A + n_mismatch_C + n_mismatch_G + n_mismatch_T;

        // check if this position has enough depth
        if(N < MIN_DEPTH) {
            continue;
        }
        ++L; // increment number of positions

        // compute D_l for this locus l
        N_Nminus1 = N*(N-1);
        D_l = (N_Nminus1 - (n_match*(n_match-1)) - (n_mismatch_A*(n_mismatch_A-1)) - (n_mismatch_C*(n_mismatch_C-1)) - (n_mismatch_G*(n_mismatch_G-1)) - (n_mismatch_T*(n_mismatch_T-1)))/N_Nminus1;
        cout << chrom << '\t' << pos << '\t' << D_l << endl;
        pi += D_l;
    }


    // output pi
    cout << "L then PI\t" << L << '\t';
    if(L == 0) {
        cout << "0";
    } else {
        pi /= L;
        cout << pi;
    }
    cout << endl;
}
