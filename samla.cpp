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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
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

static string       references_file;
static string       output_file;
static bool         opt_stdio = false;
static int32_t      opt_debug = 1;
static int32_t      debug_progress = 100000;
static int64_t      opt_progress = 0; // 1000000;
static const char tab = '\t';
static const string sep = "\t";

static void usage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "");

ostream& operator<<(ostream& out, const vector<string>& v) {
    // ------ print string vector as [i]string[,[i+1]string ...]
    // ------ e.g., [0]firststring,[1]secondstring,[2]thirdstring
    size_t i = 0;
    for (vector<string>::const_iterator p = v.begin(); p != v.end(); ++i) {
        out << "[" << i << "]" << *p;
        if (++p != v.end()) out << ",";
    }
    return(out);
}

class VcfStripmine {
    // Manage a collection of VCF files to "strip" entries off of them in order.
    // Return a vector of vcflib::Variant entries for each position we encounter.
    // If they are not in the same order, strange things will happen.
    // TODO:
    // - list of references passed in, otherwise implied by first VCF
    // - handle SNPs
    // - handle indels
    // - check for sorted order
    // - return a window's worth of variants
    // - handle direct access
    
    // TODO: friend ostream& operator<<(ostream& out, VcfStripmine& vsm);
    public:
        typedef vector<Variant const * > VariantConstPtr_vector;
        // ------------------------ Helper classes
        class VcfStrip {
            public:
                // Switched to dynamically allocated members to avoid creating an implied copy
                // constructor when creating the member in the map, maybe I could clean this up as in 
                // http://stackoverflow.com/questions/15948695/stl-map-insertion-copy-constructor
                VariantCallFile * vcf;
                Variant         * var;
                //Variant         * last_var;  // TODO: do I need this?
                string            vcf_filename;
                bool              is_open;
                bool              is_allocated;
            VcfStrip(void) : is_open(false) { }
            ~VcfStrip(void) {
                if (is_allocated) {
                    delete vcf;
                    delete var;
                    //delete last_var;
                }
            }
            bool get_next_variant(void) {
                // Call vcflib's getNextVariant() within this wrapper.  Currently it:
                // 1. skips indels
                // 2. sets is_open = false for the strip if no more variants
                while (vcf->getNextVariant(*var)) {
                    if (VcfStripmine::is_indel(*var)) {
                        if (DEBUG(1)) cerr << ".get_next_variant: skipped indel at " << vcf_filename << " : " << var->sequenceName << " : " << var->position << endl;
                        continue;
                    }
                    return(true);
                }
                if (DEBUG(1)) cerr << ".get_next_variant: no more variants in " << vcf_filename << endl;
                is_open = false;
                return(false);
            }
        };
        typedef map<string, VcfStrip>     Stripmine;
        typedef Stripmine::iterator       StripmineI;
        typedef Stripmine::const_iterator StripmineCI;
        // ------------------------ Public variables
        string    first_vcf; // first vcf file provided to the class
        string    first_ref; // first reference, inferred or provided to initiate()
        string    ref;       // current reference
        size_t    ref_index; // index of current reference in the_refs[];
        long      pos;       // current position
        string    last_ref;  // last reference
        long      last_pos;  // last position

        vector<string> the_refs;  // the collection of reference sequences; TODO see load_references() 

        Stripmine      the_mine;  // the collection of VCF files and their variants that we are stripmining

        // ------------------------ Xtor, Dtor

        VcfStripmine(void) : first_vcf(""), first_ref(""), ref(""), ref_index(0), pos(0), last_ref(""), last_pos(0) { }
        ~VcfStripmine(void) { }

        // ------------------------ Const methods

        int size(void) const { return(the_mine.size()); }
        int open_size(void) const { 
            int ans = 0;
            for (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) if (p->second.is_open) ++ans;
            return ans;
        }
        bool is_open(void) const { return(size() > 0); }
        void print(bool details = true) const {
            cerr << "VcfStripmine.print:" << endl << tab << "first_vcf:" << first_vcf << "  first_ref:" << first_ref 
                << " ref:" << ref << "  pos:" << pos << " last_ref:" << last_ref << "  last_pos:" << last_pos << endl;
            if (details) show_pos(".print", true);
        }
        void show_pos(const string& msg = "", bool details = false) const {
            cerr << "VcfStripmine.show_pos: " << msg << endl;
            cerr << tab << "VCF" << tab << "REF" << tab << "POS" << endl;
            cerr << tab << "*class*" << tab << ref << tab << pos << tab << "last_ref:" << last_ref << tab << "last_pos:" << last_pos << endl;
            for (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open)
                    cerr << tab << p->first << tab << "is_open = false" << endl;
                else {
                    cerr << tab << p->first << tab << p->second.var->sequenceName << tab << p->second.var->position << endl;
                    if (details) cerr << *p->second.var << endl;
                }
            }
        }
        bool error_false(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "") const {
            cerr << "VcfStripmine" << msg << msg2 << msg3 << msg4 << endl;
            return false;
        }
        bool error_exit(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "") const {
            cerr << "VcfStripmine" << msg << msg2 << msg3 << msg4 << endl;
            exit(EXIT_FAILURE);  // so not *technically* a const method
        }

        // ------------------------ Non-const Methods

        // TODO: allow for loading references from SAM-format lines 
        // TODO: allow for loading references from a BAM file
        // TODO: allow for references to be implied by the first VCF
        bool load_references(const string& filename) {
            ifstream refs(filename.c_str());
            if (! refs.is_open()) return(error_false(".load_references: Couldn't open reference file ", filename));
            string s;
            while (getline(refs, s)) {
                the_refs.push_back(s);
            }
            if (DEBUG(1)) cerr << "VcvStripmine.load_references: " << the_refs.size() << " reference sequences loaded from " << filename << endl;
            if (! the_refs.size()) return(error_false(".load_references: Loaded **0** reference sequence names from ", filename));
            return true;
        }
        bool add(char * filename) {
            string s(filename);
            return add(s);
        }
        bool add(string& filename) {  // (const string& filename) {
            if (the_mine.find(filename) != the_mine.end()) return(error_false(".add: VCF file already in the mine: ", filename));
            // set up a strip for this VCF file
            the_mine[filename].vcf = new VariantCallFile;
            the_mine[filename].vcf->open(filename);
            if (! the_mine[filename].vcf->is_open()) return(error_false(".add: couldn't open VCF file: ", filename));
            the_mine[filename].vcf_filename = filename;
            the_mine[filename].is_open = true;
            the_mine[filename].var = new Variant(*the_mine[filename].vcf);
            the_mine[filename].is_allocated = true;

            // or: the_mine[filename].var = new Variant; the_mine[filename].var->setVariantCallFile(the_mine[filename].vcf);

            // do we need last_var?  the_mine[filename].last_var = new Variant(*the_mine[filename].vcf);

            if (! first_vcf.size()) first_vcf = filename;
            return true;
        }
        bool remove(const string& filename) {
            StripmineI p = the_mine.find(filename);
            if (p == the_mine.end()) return false;
            the_mine.erase(p);
            if (filename == first_vcf) first_vcf = "";
            return true;
        }
        bool initiate(void) {
            // TODO:  check for already-initiated VCFs
            // TODO:  handle reinitiation?
            if (the_mine.size() == 0) return(error_false(".initiate: no vcf files add()ed"));
            // process the first variant in each VCF
            // TODO: someday: can combine the following and this loop if I change the_refs[] a bit
            int n = 0;
            for (StripmineI p = the_mine.begin(); p != the_mine.end(); ++p)
                if (p->second.get_next_variant()) ++n;
            if (! n) return(error_false(".initiate: couldn't find variants in any of the files in the mine"));
            // ref_index = 0 should be a no-op, it is initialised to 0 by our ctor
            for (ref_index = 0; ref_index < the_refs.size(); ++ref_index) {
                for (StripmineI p = the_mine.begin(); p != the_mine.end(); ++p) {
                    if (the_refs[ref_index] == p->second.var->sequenceName) {
                        ref = the_refs[ref_index];
                        break;
                    }
                }
                if (! ref.empty()) break;
            }
            if (ref.empty()) return(error_false(".initiate: couldn't find variants for any of the supplied references"));
            pos = 0;
            assert(! ref.empty() && pos == 0);
            if (DEBUG(1)) show_pos(".initiate, ready to return");
            return true;
        }
    private:
        bool find_next_pos(void) {
            // we think there might be at least one variant left on this ref, and all variants are unharvested
            // return true if we found one, return false if we didn't (false means we will have to look at the next ref)
            long new_pos = pos;
            for (StripmineI p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open) continue;
                if (p->second.var->sequenceName == ref && p->second.var->position > pos) {
                    if (new_pos == pos || ((p->second.var->position - pos) < (new_pos - pos)))
                        new_pos = p->second.var->position;
                }
            }
            if (new_pos == pos) {
                if (DEBUG(1)) cerr << ".find_next_pos: could not find a new variant" << endl;
                return(false);
            }
            last_ref = ref;
            last_pos = pos;
            pos = new_pos;
            return(true);
        }
        bool find_next_pos_new_ref(void) {
            // we don't have a ref or we have to switch to a new one, and all variants are unharvested
            if (DEBUG(2)) print(true);
            string new_ref;
            long new_pos = 0;
            for (++ref_index; ref_index < the_refs.size(); ++ref_index) {
                for (StripmineI p = the_mine.begin(); p != the_mine.end(); ++p) {
                    if (! p->second.is_open) continue;
                    if (the_refs[ref_index] == p->second.var->sequenceName) {
                        new_ref = the_refs[ref_index];
                        if (new_pos == 0 || new_pos > p->second.var->position) new_pos = p->second.var->position;
                        break;
                    }
                }
                if (! new_ref.empty()) break;
            }
            if (new_ref.empty() || new_pos == 0) {
                if (DEBUG(1)) cerr << ".find_next_pos_new_ref: could not find a new variant" << endl;
                if (DEBUG(2)) print(true);
                return(false);
            }
            last_ref = ref;
            last_pos = pos;
            ref = new_ref;
            pos = new_pos;
            if (DEBUG(2)) print(true);
            return(true);
        }
        void refresh_variants(void) {
            StripmineI p;
            // called to make all variants unharvested; we assume that any variants
            // that match the current ref and pos were harvested in the last round, so get
            // the next variant in each such vcf, or close it
            for (p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open) continue;
                if (p->second.var->sequenceName == ref && p->second.var->position == pos) {
                    if (DEBUG(2)) cerr << ".refresh_variants: getting new variant from " << p->first << endl;
                    p->second.get_next_variant();
                }
            }
        }
        bool get_at(VariantConstPtr_vector& vars, const string& this_ref, const long this_pos) {
            // get all variants for ref rf, position ps; currently only works if there are already VCFs positioned there
            int i = 0;
            for (StripmineCI p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (p->second.var->sequenceName == this_ref && p->second.var->position == this_pos) {
                    vars.push_back(p->second.var);
                    ++i;
                }
            }
            if (i == 0) {
                cerr << ".get_at: failed to find variants at " << this_ref << " : " << this_pos << endl;
                return(false);
            }
            return(true);
        }
    public:
        bool get(VariantConstPtr_vector& vars) {
            // Returns the next set of variants in the mine.  Upon exit, vars[] holds the (immutable)
            // variants, ref and pos of the class point to the location of the variants being returned,
            // and last_ref and last_pos hold the location of the last variant returned.
            //
            // On entry, the_mine[].var holds either the last variant (.sequenceName and .position match ref and pos) 
            // or one that may be upcoming.  If pos == 0, then the class was just initiated, ref has been set to the
            // reference for the first variant (checked for reference name only), and all the_mine[].var 
            // variants are unharvested
            if (ref.empty()) error_exit("VcfStripmine.get: must call .initiate() before first call to .get()");
            vars.clear();
            if (pos == 0) { 
                // just initialised, all variants unharvested, at least one variant holds a variant in the reference named in ref
                if (! find_next_pos())
                    return(false);
            } else {  
                // need to find a new variant, either on this ref or on the next ref holding a variant
                refresh_variants();
                if (! find_next_pos() && ! find_next_pos_new_ref()) {
                    return(false);
                }
            }
            // at this point, ref and pos have been successfully set to point to the location of the next variant(s)
            get_at(vars, ref, pos);
            if (DEBUG(1)) cerr << ".get: returning " << vars.size() << " variants at " << ref << " : " << pos << endl;
            assert(vars.size());
            return(true);
        }

        // ------------------------ "Helper" methods... nothing class-specific, only here to restrict their scope

        static bool is_indel(Variant& var) {
            // I think the only sure way to tell is if the ref and alt alleles are of different lengths
            size_t diff_lengths = 0;
            for (size_t i = 1; i < var.alleles.size(); ++i) 
                if (var.alleles[0].length() != var.alleles[i].length()) ++diff_lengths;
            return(diff_lengths > 0);
        }

    // end of public: for Class VcfStripmine

}; // end of Class VcfStripmine



//-------------------------------------

int 
main(int argc, char* argv[])
{
    VcfStripmine vcfmine;

    //----------------- Command-line options

    enum { o_references, o_output, o_stdio, o_debug, o_progress, o_help };

    CSimpleOpt::SOption smorgas_options[] = {
        { o_references, "-r",           SO_REQ_SEP },
        { o_references, "--references", SO_REQ_SEP },
        { o_output,     "-o",           SO_REQ_SEP },
        { o_output,     "--output",     SO_REQ_SEP },
        { o_stdio,      "-",            SO_NONE }, 
        { o_debug,      "--debug",      SO_REQ_SEP },
        { o_progress,   "--progress",   SO_REQ_SEP },
        { o_help,       "--help",       SO_NONE },
        { o_help,       "-?",          SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, smorgas_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            //cerr << NAME + " invalid argument '" << args.OptionText() << "'" << endl;
            usage("invalid argument '", args.OptionText(), "'");
        }
        switch (args.OptionId()) {
            case o_help:     
                usage(); break;
            case o_references:   
                references_file = args.OptionArg(); break;
            case o_output:   
                output_file = args.OptionArg(); break;
            case o_stdio:    
                opt_stdio = true; break;
            case o_debug:    
                opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug; break;
            case o_progress: 
                opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress; break;
            default:         
                usage("invalid option '", args.OptionText(), "'"); break;
        }
    }

    if (args.FileCount() == 0) usage("At least one VCF file must be specified as input");

    // read in reference sequence names
    if (references_file.empty()) usage("Give file holding reference seq names in order with -r/--references");
    vcfmine.load_references(references_file);

    // add VCF file to the mine
    for (int i = 0; i < args.FileCount(); ++i) {
        vcfmine.add(args.File(i));
    }

    if (output_file.empty()) {
        output_file = "/dev/stdout"; 
    } else {
        usage("********** output file not yet supported, sorry");
    }

    // other options to initiate
    if (DEBUG(1) and ! opt_progress) opt_progress = debug_progress;

    // initiate the stripmine, which fetches the first variant from each file
    if (! vcfmine.initiate()) usage("Could not initiate stripmine");

    //vcfmine.print();

    VcfStripmine::VariantConstPtr_vector vars;

    while (vcfmine.get(vars)) {

        for (VcfStripmine::VariantConstPtr_vector::const_iterator p = vars.begin(); p != vars.end(); ++p) {
        }

    }

    return(EXIT_SUCCESS);
}


//-----------------


static void
//usage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "")
usage(const string& msg, const string& msg2, const string& msg3, const string& msg4)
{
    if (! msg.empty()) cerr << msg << msg2 << msg3 << msg4 << endl;
    cerr << endl;
    cerr << "Usage:   " << NAME << " [options] -r refnames.txt <in1.vcf> [ in2.vcf ... ]" << endl;
    cerr << endl;
    cerr << "Collect VCF results and produce a consensus list of variants." << endl;
    cerr << endl;
    cerr << "NOTE: This command is under active development." << endl;
    cerr << endl;
    cerr << "     -r FILE | --references FILE  file containing reference names in order [REQUIRED]" << endl;
    cerr << "                                  should be in the same order as the VCF files" << endl;
    cerr << "     -o FILE | --output FILE      output file name [default is stdout]" << endl;
    cerr << "     --debug INT                  debug info level INT [" << opt_debug << "]" << endl;
    cerr << "     --progress INT               print variants processed mod INT [" << opt_progress << "]" << endl;
    cerr << endl;
    cerr << "     -? | --help                  help" << endl;
    cerr << endl;
    exit(EXIT_FAILURE);
}


