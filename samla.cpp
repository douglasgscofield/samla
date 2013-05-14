// samla.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// samla [Swedish]: (v) to collect
// Collect VCF files containing variants from multiple sources and interpret
// the results
//
// Uses the vcflib (https://github.com/ekg/vcflib) library for reading VCF
// files.
//

// CHANGELOG
//
//

// TODO
// --- open and read VCF file
// --- " " " multiple VCF files
// --- output messages combining info from variants across VCF files
// --- deal with indels
// --- combine likelihoods
//

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>
#define IF_DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define DEBUG(__lvl__) (opt_debug >= __lvl__)

// vcflib
//
#include "Variant.h"

// very nice argument handling
//
#include "SimpleOpt.h"

// Std C/C++ includes, many already included by Variant.h

//#include <cstdlib>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include <vector>
//#include <string>
//#include <memory>
//#include <limits>
//#include <algorithm>
//#include <ctype.h>
//#include <stdint.h>

#define NAME "samla"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

using namespace std;
using namespace vcf;

static string       output_file;
static bool         opt_stdio = false;
static bool         opt_mappingquality = false;
static bool         opt_profile = false;
static int32_t      opt_debug = 1;
static int32_t      debug_progress = 100000;
static int64_t      opt_progress = 0; // 1000000;
static const char tab = '\t';
static const string sep = "\t";



//-------------------------------------

int 
main(int argc, char* argv[])
{
    VcfStripmine vcfmine;

    //----------------- Command-line options

    enum { o_output, o_stdio, o_debug, o_progress, o_help };

    CSimpleOpt::SOption smorgas_options[] = {
        { o_output,   "-o",          SO_REQ_SEP },
        { o_output,   "--output",    SO_REQ_SEP },
        { o_stdio,    "-",           SO_NONE }, 
        { o_debug,    "--debug",     SO_REQ_SEP },
        { o_progress, "--progress",  SO_REQ_SEP },
        { o_help,     "--help",      SO_NONE }, { o_help,     "-?",          SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, smorgas_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            cerr << NAME << " invalid argument '" << args.OptionText() << "'" << endl;
            return usage();
        }
        switch (args.OptionID()) {
            case o_help:     return usage(true); break;
            case o_output:   output_file = args.OptionArg(); break;
            case o_stdio:    opt_stdio = true; break;
            case o_debug:    opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug; break;
            case o_progress: opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress; break;
            default:         cerr << NAME << " invalid option '" << args.OptionText() << "'" << endl; return usage(true); break;
        }
    }
    if (args.FileCount() == 0) {
        cerr << NAME << " requires at least one VCF file specified as input" << endl;
        return usage();
    }

    if (DEBUG(1) and ! opt_progress) opt_progress = debug_progress;

    for (int i = 0; i < args.FileCount(); ++i) {
        vcfmine.add(args.File(i));
    }

    if (output_file.empty()) {
        output_file = "/dev/stdout"; 
    } else {
        cerr << NAME << ": output file not yet supported" << endl;
        return EXIT_FAILURE;
    }

    if (! vcfmine.initiate()) return EXIT_FAILURE;

    vector<Variant> vars = vcfmine.get();

    while (vars.size()) {

        vars = vcfmine.get();
    }
}


class VcfStripmine {
    // Manage a collection of VCF files to "strip" entries off of them in order.
    // Return a vector of vcflib::Variant entries for each position we encounter.
    // If they are not in the same order, strange things will happen.
    // TODO:
    // - handle SNPs
    // - handle indels
    // - check for sorted order
    // - return a window's worth of variants
    // - handle direct access
    friend std::ostream& operator<<(std::ostream& out, VcfStripmine& vsm);
    public:
        // ------------------------ Helper classes
        class VcfStrip {
            public:
                VariantCallFile vcf;
                Variant         var;
                Variant         last_var;
                bool            is_open;
            VcfStrip(void) : is_open(false) { }
        };
        typedef std::map<std::string, VcfStrip>  Stripmine;
        typedef Stripmine::iterator              StripmineI;
        typedef Stripmine::const_iterator        StripmineCI;
        // ------------------------ Public variables
        string    reference; // current reference
        long      pos;     // current position
        string    last_reference; // last reference
        long      last_pos;     // last position
        Stripmine the_mine;
        // ------------------------ Xtor, Dtor
        VcfStripmine(void) : num_files(0), reference(""), pos(0), last_reference(""), last_pos(0) { }
        ~VcfStripmine(void) { }
        // ------------------------ Const methods
        int size(void) const { return(the_mine.size()); }
        int open_size(void) const { 
            int ans = 0;
            for (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) if (p->second.is_open) ++ans;
            return ans;
        }
        bool is_open(void) const { return(size() > 0); }
        void print(void) const {
            for (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) {
                cerr << "VcfStripmine.print:" << p->first << ":vcf  : " << p->second.vcf << endl;
                cerr << "VcfStripmine.print:" << p->first << ":var  : " << p->second.var << endl;
                cerr << "VcfStripmine.print:" << p->first << ":last_: " << p->second.last_var << endl;
                cerr << "VcfStripmine.print:" << p->first << ":is_op: " << p->second.is_open << endl;
            }
        }
        void show_pos(const std::string& msg = "") {
            while (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open) 
                    cerr << "VcfStripmine" << msg << ": " << p->first << " is_open = false" << endl;
                else 
                    cerr << "VcfStripmine" << msg << ": " << p->first << " @ " << p->second.var.sequenceName << " : " 
                        << p->second.var.position << endl;
            }
        }
        // ------------------------ Methods
        bool add(const std::string& filename) {
            if (the_mine.find(filename) != the_mine.end()) {
                cerr << "VcfStripmine.add: VCF file already in the mine: " << filename << endl;
                return false;
            }
            VariantCallFile vcffile;
            vcffile.open(filename);
            if (! vcffile.is_open()) { 
                cerr << "VcfStripmine.add: couldn't open VCF file: " << filename << endl;
                return false;
            }
            the_mine[filename].vcf = vcffile;
            the_mine[filename].var(the_mine[filename].vcf);
            the_mine[filename].last_var(the_mine[filename].vcf);
            the_mine[filename].is_open = true;
            return true;
        }
        bool remove(const std::string& filename) {
            StripmineI p = the_mine.find(filename);
            if (p == the_mine.end()) return false;
            the_mine.erase(p);
            return true;
        }
        bool initiate(void) {
            // TODO:  check for already-open VCFs
            // TODO:  handle reinitiation
            if (the_mine.size() == 0) {
                cerr << "VcfStripmine.initiate: no files in the mine" << endl;
                return false;
            }
            int n_init = 0;
            while (StripmineI p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (p->second.vcf.getNextVariant(p->second.var)) ++n_init;
                else p->second.is_open = false;
            }
            show_pos(".initiate");
            if (! n_init) {
                cerr << "VcfStripmine.initiate: couldn't find variants in any of the files in the mine" << endl;
                return false;
            }
            return true;
        }
        std::vector<Variant> get() {
            std::vector<Variant> ans;
            if (reference == "" || pos == 0) {
                cerr << "VcfStripmine.get: must call .initiate() before first call to .get()" << endl;
                exit(EXIT_FAILURE);
            }
            if (DEBUG(1)) show_pos(".get start");
            std::string
            if (DEBUG(1)) show_pos(".get end");
        }
}; // Class VcfStripmine

//-----------------


static int
usage(bool longer = false)
{
    cerr << endl;
    cerr << "Usage:   " << NAME << " [options] <in1.vcf> [ in2.vcf ... ]" << endl;
    cerr << endl;
    cerr << "Collect VCF results and produce a consensus list of variants." << endl;
    cerr << endl;
    cerr << "NOTE: This command is very much a work in progress." << endl;
    cerr << endl;
    cerr << "     -o FILE | --output FILE   output file name [default is stdout]" << endl;
    cerr << "     -? | --help               output file name [default is stdout]" << endl;
    cerr << "     --debug INT               debug info level INT [" << opt_debug << "]" << endl;
    cerr << "     --progress INT   print reads processed mod INT [" << opt_progress << "]" << endl;
    cerr << endl;
    return EXIT_FAILURE;
}


