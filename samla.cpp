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
// -x- open and read VCF file
// -x- " " " multiple VCF files
// --- output messages combining info from variants across VCF files
// --- deal with indels
// --- combine likelihoods
// --- return Variants in order of VCFs
//

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>
#define IF_DEBUG(__lvl__) if (opt_debug >= __lvl__)
#define DEBUG(__lvl__) (opt_debug >= __lvl__)

// vcflib
//
#include "src/Variant.h"

// very nice argument handling
//
#include "SimpleOpt.h"

// Std C/C++ includes, many already included by Variant.h

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <deque>
#include <vector>
#include <string>
//#include <memory>
//#include <limits>
#include <algorithm>
//#include <ctype.h>
//#include <stdint.h>

#define NAME "samla"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

using namespace std;
using namespace vcf;

static string       references_file;
static string       output_file;
static string       samla_method;
static string       vcf_genomic;
static string       vcf_wga;
static string       vcf_all;
static const size_t max_gwa_window_size = 50;
static size_t       opt_gwa_window_size = 5;
static double       opt_gwa_quality = 30.0;
static bool         opt_stdio = false;
static int32_t      opt_debug = 1;
static int32_t      debug_progress = 100000;
static int64_t      opt_progress = 0; // 1000000;
static const char tab = '\t';
static const string sep = "\t";

static void exitmessage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "", const string& msg5 = "");
static void exitusage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "", const string& msg5 = "");


template<class T>
ostream& operator<<(ostream& out, const vector<T>& v) {
    // ------ print vector as [i]val[,[i+1]val ...]
    // ------ e.g., [0]firstval,[1]secondval,[2]thirdval
    size_t i = 0;
    for (typename vector<T>::const_iterator p = v.begin(); p != v.end(); ++i) {
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
    //
    // TODO:
    // - what about multiple samples per vcf line?
    
    // TODO: friend ostream& operator<<(ostream& out, VcfStripmine& vsm);

    public:

        typedef vector<Variant const * >         VariantConstPtr_vector;
        typedef vector<Variant * >               VariantPtr_vector;
        typedef vector<VariantCallFile const * > VariantCallFileConstPtr_vector;
        typedef vector<VariantCallFile * >       VariantCallFilePtr_vector;


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
                            if (DEBUG(1)) 
                                cerr << ".get_next_variant(): skipped indel at " << vcf_filename 
                                    << " : " << var->sequenceName << " : " << var->position << endl;
                            continue;
                        }
                        return(true);
                    }
                    if (DEBUG(1)) 
                        cerr << ".get_next_variant(): no more variants in " << vcf_filename << endl;
                    is_open = false;
                    return(false);
                }

        };

        typedef map<string, VcfStrip>     Stripmine;  // maps VCF filenames to VcfStrip instances


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


        VcfStripmine(void) : first_vcf(""), first_ref(""), ref(""), ref_index(0), pos(0), 
                             last_ref(""), last_pos(0) { }
        ~VcfStripmine(void) { }


        // ------------------------ Const methods


        int size(void) const { 
            return(the_mine.size()); 
        }

        int open_size(void) const { 
            int ans = 0;
            for (Stripmine::const_iterator p = the_mine.begin(); p != the_mine.end(); ++p) 
                if (p->second.is_open) 
                    ++ans;
            return ans;
        }

        bool is_open(void) const { 
            return(size() > 0); 
        }

        void print(bool details = true) const {
            cerr << ".print():" << endl << tab << "first_vcf:" << first_vcf 
                << "  first_ref:" << first_ref << " ref:" << ref << "  pos:" << pos 
                << " last_ref:" << last_ref << "  last_pos:" << last_pos << endl;
            if (details) show_pos(".print", true);
        }

        void show_pos(const string& msg = "", bool details = false) const {
            cerr << ".show_pos(): " << msg << endl;
            fprintf(stderr, "\t%-30s\t%-20s\t%-10s\n", "VCF", "REF", "POS");
            fprintf(stderr, "\t%-30s\t%-20s\t%10ld\tlast_ref:%s\tlast_pos:%ld\n", "*class*", 
                    ref.c_str(), pos, last_ref.c_str(), last_pos);
            for (Stripmine::const_iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open)
                    fprintf(stderr, "\t%-30s\t%-20s\n", p->first.c_str(), "is_open = false");
                else {
                    fprintf(stderr, "\t%-30s\t%-20s\t%10ld\n", p->first.c_str(), 
                            p->second.var->sequenceName.c_str(), p->second.var->position);
                    if (details) 
                        cerr << *p->second.var << endl;
                }
            }
        }

        bool false_error(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "") const {
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
            if (! refs.is_open()) 
                return(false_error(".load_references(): Couldn't open reference file ", filename));
            string s;
            while (getline(refs, s))
                the_refs.push_back(s);
            if (DEBUG(1)) 
                cerr << ".load_references(): " << the_refs.size() 
                    << " reference sequences loaded from " << filename << endl;
            if (! the_refs.size()) 
                return(false_error(".load_references(): Loaded **0** reference sequence names from ", filename));
            return true;
        }

        bool add(char * filename) {
            string s(filename);
            return add(s);
        }

        bool add(string& filename) {  // (const string& filename) {
            if (the_mine.find(filename) != the_mine.end()) 
                return(false_error(".add(): VCF file already in the mine: ", filename));
            // set up a strip for this VCF file
            the_mine[filename].vcf = new VariantCallFile;
            the_mine[filename].vcf->open(filename);
            if (! the_mine[filename].vcf->is_open()) 
                return(false_error(".add(): couldn't open VCF file: ", filename));
            the_mine[filename].vcf_filename = filename;
            the_mine[filename].is_open = true;
            the_mine[filename].var = new Variant(*the_mine[filename].vcf);
            the_mine[filename].is_allocated = true;

            // or: the_mine[filename].var = new Variant; the_mine[filename].var->setVariantCallFile(the_mine[filename].vcf);

            // do we need last_var?  the_mine[filename].last_var = new Variant(*the_mine[filename].vcf);

            if (! first_vcf.size()) 
                first_vcf = filename;
            return true;
        }

        bool remove(const string& filename) {
            Stripmine::iterator p = the_mine.find(filename);
            if (p == the_mine.end()) 
                return false;
            the_mine.erase(p);
            if (filename == first_vcf) 
                first_vcf = "";
            return true;
        }

        bool initiate(void) {
            // TODO:  check for already-initiated VCFs
            // TODO:  handle reinitiation?
            if (the_mine.size() == 0) 
                return(false_error(".initiate(): no vcf files add()ed"));
            // process the first variant in each VCF
            // TODO: someday: can combine the following and this loop if I change the_refs[] a bit
            int n = 0;
            for (Stripmine::iterator p = the_mine.begin(); p != the_mine.end(); ++p)
                if (p->second.get_next_variant()) 
                    ++n;
            if (! n) 
                return(false_error(".initiate(): couldn't find variants in any of the files in the mine"));
            // ref_index = 0 should be a no-op, it is initialised to 0 by our ctor
            for (ref_index = 0; ref_index < the_refs.size(); ++ref_index) {
                for (Stripmine::iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                    if (the_refs[ref_index] == p->second.var->sequenceName) {
                        ref = the_refs[ref_index];
                        break;
                    }
                }
                if (! ref.empty()) 
                    break;
            }
            if (ref.empty()) 
                return(false_error(".initiate(): couldn't find variants for any of the supplied references"));
            pos = 0;
            assert(! ref.empty() && pos == 0);
            if (DEBUG(3)) 
                show_pos(".initiate(), ready to return");
            return true;
        }

    private:

        bool find_next_pos(void) {
            // we think there might be at least one variant left on this ref, and all variants are unharvested
            // return true if we found one, return false if we didn't (false means we will have to look at the next ref)
            long new_pos = pos;
            for (Stripmine::iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open) 
                    continue;
                if (p->second.var->sequenceName == ref && p->second.var->position > pos) {
                    if (new_pos == pos || ((p->second.var->position - pos) < (new_pos - pos)))
                        new_pos = p->second.var->position;
                }
            }
            if (new_pos == pos) {
                if (DEBUG(1)) 
                    cerr << ".find_next_pos(): could not find a new variant" << endl;
                return(false);
            }
            last_ref = ref;
            last_pos = pos;
            pos = new_pos;
            return(true);
        }

        bool find_next_pos_new_ref(void) {
            // we don't have a ref or we have to switch to a new one, and all variants are unharvested
            if (DEBUG(3)) 
                print(true);
            string new_ref;
            long new_pos = 0;
            for (++ref_index; ref_index < the_refs.size(); ++ref_index) {
                for (Stripmine::iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                    if (! p->second.is_open) 
                        continue;
                    if (the_refs[ref_index] == p->second.var->sequenceName) {
                        new_ref = the_refs[ref_index];
                        if (new_pos == 0 || new_pos > p->second.var->position) 
                            new_pos = p->second.var->position;
                        break;
                    }
                }
                if (! new_ref.empty()) 
                    break;
            }
            if (new_ref.empty() || new_pos == 0) {
                if (DEBUG(1)) 
                    cerr << ".find_next_pos_new_ref(): could not find a new variant" << endl;
                if (DEBUG(3)) 
                    print(true);
                return(false);
            }
            last_ref = ref;
            last_pos = pos;
            ref = new_ref;
            pos = new_pos;
            if (DEBUG(3)) 
                print(true);
            return(true);
        }

        void refresh_variants(void) {
            Stripmine::iterator p;
            // called to make all variants unharvested; we assume that any variants
            // that match the current ref and pos were harvested in the last round, so get
            // the next variant in each such vcf, or close it
            for (p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (! p->second.is_open) 
                    continue;
                if (p->second.var->sequenceName == ref && p->second.var->position == pos) {
                    if (DEBUG(3)) 
                        cerr << ".refresh_variants: getting new variant from " << p->first << endl;
                    p->second.get_next_variant();
                }
            }
        }

        bool get_at(VariantConstPtr_vector& vars, const string& this_ref, const long this_pos) {
            // get all variants for ref rf, position ps; currently only works if there are already VCFs positioned there
            int i = 0;
            for (Stripmine::const_iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                if (p->second.var->sequenceName == this_ref && p->second.var->position == this_pos) {
                    vars.push_back(p->second.var);
                    ++i;
                }
            }
            if (i == 0) {
                cerr << ".get_at: failed to find Variant at " << this_ref << " : " << this_pos << endl;
                return(false);
            }
            return(true);
        }

    public:

        bool get_files(VariantCallFileConstPtr_vector& f) {
            if (the_mine.size() == 0) 
                return(false_error(".get_files: no vcf files add()ed"));
            f.clear();
            for (Stripmine::const_iterator p = the_mine.begin(); p != the_mine.end(); ++p) {
                f.push_back(p->second.vcf);
            }
            return(true);
        }
        bool get(VariantConstPtr_vector& vars) {
            // Returns the next set of variants in the mine.  Upon exit, vars[] holds the (immutable)
            // variants, ref and pos of the class point to the location of the variants being returned,
            // and last_ref and last_pos hold the location of the last variant returned.
            //
            // On entry, the_mine[].var holds either the last variant (.sequenceName and .position match ref and pos) 
            // or one that may be upcoming.  If pos == 0, then the class was just initiated, ref has been set to the
            // reference for the first variant (checked for reference name only), and all the_mine[].var 
            // variants are unharvested
            if (ref.empty()) 
                error_exit("VcfStripmine.get: must call .initiate() before first call to .get()");
            vars.clear();
            if (pos == 0) { 
                // 1+ the_mine[].var hold a variant in reference ref, so returning false here is unlikely
                if (! find_next_pos())
                    return(false);
            } else {  
                refresh_variants();
                if (! find_next_pos() && ! find_next_pos_new_ref()) {  // must be called in this order
                    return(false);
                }
            }
            // ref and pos point to the location of the next variant(s); perhaps we should unroll stacked
            // alternative alleles for each variant? VCF is a pretty flexible format...
            get_at(vars, ref, pos);
            if (DEBUG(3)) {
                cerr << ".get(): returning " << vars.size() << " Variant at pos " 
                    << ref << " : " << pos << " : N-samples";
                for (size_t i = 0; i < vars.size(); ++i)
                    cerr << " " << vars[i]->getNumSamples();
                cerr << endl;
            }
            assert(vars.size());
            return(true);
        }

        // ------------------------ "Helper" methods... nothing class-specific, only here to restrict their scope

        static bool is_indel(Variant& var) {
            // I think the only sure way to tell is if the ref and alt alleles are of different lengths
            size_t diff_lengths = 0;
            for (size_t i = 1; i < var.alleles.size(); ++i) 
                if (var.alleles[0].length() != var.alleles[i].length()) 
                    ++diff_lengths;
            return(diff_lengths > 0);
        }

    // end of public: for Class VcfStripmine

}; // end of Class VcfStripmine


//-------------------------------------


void
processCommandLine(int argc, char* argv[], VcfStripmine& stripmine) {

    enum { o_references, o_method, 
           o_vcf_genomic, o_vcf_wga, o_vcf_all, o_gwa_window, o_gwa_quality,
           o_output, o_stdio, o_debug, o_progress, o_help };

    CSimpleOpt::SOption smorgas_options[] = {
        { o_references,  "-r",            SO_REQ_SEP },  // file of reference sequence names
        { o_references,  "--references",  SO_REQ_SEP },
        { o_method,      "-m",            SO_REQ_SEP },  // method for variant combining
        { o_method,      "--method",      SO_REQ_SEP },

        { o_vcf_genomic, "--vcf-genomic", SO_REQ_SEP },
        { o_vcf_genomic, "-g",            SO_REQ_SEP },
        { o_vcf_wga,     "--vcf-wga",     SO_REQ_SEP },
        { o_vcf_wga,     "-w",            SO_REQ_SEP },
        { o_vcf_all,     "--vcf-all",     SO_REQ_SEP },
        { o_vcf_all,     "-a",            SO_REQ_SEP },
        { o_gwa_window,  "--gwa-window",  SO_REQ_SEP },
        { o_gwa_quality, "--gwa-quality",  SO_REQ_SEP },

        { o_output,      "-o",            SO_REQ_SEP },  // output file
        { o_output,      "--output",      SO_REQ_SEP },
        { o_stdio,       "-",             SO_NONE }, 
        { o_debug,       "--debug",       SO_REQ_SEP },
        { o_progress,    "--progress",    SO_REQ_SEP },
        { o_help,        "--help",        SO_NONE },
        { o_help,        "-?",            SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    CSimpleOpt args(argc, argv, smorgas_options);

    while (args.Next()) {
        if (args.LastError() != SO_SUCCESS) {
            //cerr << NAME + " invalid argument '" << args.OptionText() << "'" << endl;
            exitusage("invalid argument '", args.OptionText(), "'");
        }
        switch (args.OptionId()) {
            case o_help:     
                exitusage(); break;
            case o_references:   
                references_file = args.OptionArg(); break;
            case o_method:   
                samla_method = args.OptionArg(); break;
            case o_vcf_genomic:   
                vcf_genomic = args.OptionArg(); break;
            case o_vcf_wga:   
                vcf_wga = args.OptionArg(); break;
            case o_vcf_all:   
                vcf_all = args.OptionArg(); break;
            case o_gwa_window:   
                opt_gwa_window_size = args.OptionArg() ? atoi(args.OptionArg()) : opt_gwa_window_size; break;
            case o_gwa_quality:   
                opt_gwa_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_quality; break;
            case o_output:   
                output_file = args.OptionArg(); break;
            case o_stdio:    
                opt_stdio = true; break;
            case o_debug:    
                opt_debug = args.OptionArg() ? atoi(args.OptionArg()) : opt_debug; break;
            case o_progress: 
                opt_progress = args.OptionArg() ? strtoll(args.OptionArg(), NULL, 10) : opt_progress; break;
            default:         
                exitusage("invalid option '", args.OptionText(), "'"); break;
        }
    }

    // read in reference sequence names

    if (references_file.empty()) 
        exitusage("Give file holding reference seq names in order with -r/--references");
    stripmine.load_references(references_file);

    // handle methods

    if (samla_method == "gwa" || samla_method == "GWA") {

        samla_method = "gwa";
        if (opt_gwa_window_size > max_gwa_window_size)
            exitusage("--gwa-window INT too large");
        if (args.FileCount() || vcf_genomic.empty() || vcf_wga.empty() || vcf_all.empty())
            exitusage("For method 'gwa', all VCFs must be specified using --vcf-* options");
        stripmine.add(vcf_genomic);  // I assumed genomic is 0-th Variant etc. but it may not be, see TODO
        stripmine.add(vcf_wga);
        stripmine.add(vcf_all);
        if (DEBUG(1)) {
            cerr << "'gwa' method genomic VCF   : " << vcf_genomic << endl;
            cerr << "'gwa' method WGA VCF       : " << vcf_wga << endl;
            cerr << "'gwa' method all-calls VCF : " << vcf_all << endl;
        }

    } else if (samla_method == "default") {

        if (args.FileCount() == 0) 
            exitusage("At least one VCF file must be specified as input");
        // add VCF file to the mine
        for (int i = 0; i < args.FileCount(); ++i)
            stripmine.add(args.File(i));

    } else {

        exitusage("Unrecognized --method '", samla_method, "'");

    }

    if (output_file.empty())
        output_file = "/dev/stdout"; 
    else
        exitusage("********** output file not yet supported, sorry");
    //
    // other options to initiate
    if (DEBUG(1) and ! opt_progress) 
        opt_progress = debug_progress;

}


//-------------------------------------


static void
exitmessage(const string& msg, const string& msg2, const string& msg3, const string& msg4, const string& msg5) {
    if (! msg.empty()) cerr << endl << "*** " << msg << msg2 << msg3 << msg4 << msg5 << endl;
    exit(EXIT_FAILURE);
}
static void
exitusage(const string& msg, const string& msg2, const string& msg3, const string& msg4, const string& msg5) {
    if (! msg.empty()) cerr << endl << "*** " << msg << msg2 << msg3 << msg4 << msg5 << endl;
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
    cerr << "     -m METHOD | --method METHOD  use combining method 'METHOD', only 'gwa' is implemented" << endl;
    cerr << endl;
    cerr << "'gwa' method options:" <<endl;
    cerr << endl;
    cerr << "     --gwa-window INT             lookback window size for mean quality, max " << max_gwa_window_size << "[" << opt_gwa_window_size << "]" << endl;
    cerr << "     --gwa-quality FLOAT          minimum quality when combining variants [" << opt_gwa_quality << "]" << endl;
    cerr << endl;
    cerr << "For 'gwa', all VCF files must be specified using these options:" << endl;
    cerr << endl;
    cerr << "     --vcf-genomic FILE           VCF file containing genomic calls" << endl;
    cerr << "     --vcf-wga FILE               VCF file containing whole-genome-amplified calls" << endl;
    cerr << "     --vcf-all FILE               VCF file containing all (pooled) calls" << endl;
    cerr << endl;
    cerr << "     -? | --help                  help" << endl;
    cerr << endl;
    exit(EXIT_FAILURE);
}


//-------------------------------------


class LookbackWindow {  // a sliding window that holds values at the current and previous positions
    private:
        size_t        max_sz;    // maximum number of values to hold
        size_t        win_sz;    // window size of values to consider for e.g. mean (<= max_size)
        size_t        cur_sz;    // current number of values held, may be less than max_sz if just cleared.
        size_t        front_tm;  // values to ignore at the front, e.g. to exclude current position
        deque<double> buffer;    // circular buffer holding values
    public:
        LookbackWindow(size_t msz = max_gwa_window_size, size_t wsz = opt_gwa_window_size, size_t ft = 0) 
            : max_sz(msz), win_sz(wsz), cur_sz(0), front_tm(ft) { buffer.resize(max_sz); };
        size_t cur_size(void) { return (cur_sz); }
        size_t win_size(void) { return (win_sz); }
        size_t max_size(void) { return (max_sz); }
        size_t front_trim(void) { return (front_tm); }
        void   clear(void) { cur_sz = 0; buffer.clear(); buffer.resize(max_sz); }
        void   dump() { dump(win_sz, front_tm); }
        void   dump(size_t n, size_t ft) { 
            if ((n + ft) > cur_sz) n = cur_sz - ft;
            cerr << "LookbackWindow(max_sz=" << max_sz << ",win_sz=" << win_sz << ",front_tm=" << front_tm << "): ";
            cerr << "["; for (size_t i = 0; i < front_tm; ++i) { cerr << buffer[i] << ", "; } cerr << "trimmed] ";
            for (size_t i = ft; i < (n + ft); ++i) {
                cerr << buffer[i]; if (i == win_sz + ft) cerr << "|win_sz"; if (i + 1 < n + ft) cerr << ", ";
            }
            cerr << " + " << (cur_sz - n - ft) << " additional values not within " << n + ft << " of the start" << endl;
        }
        void   push(double x) {  // add a value
            if (cur_sz == max_sz) { // pop one off the end before adding
                buffer.pop_back();
                buffer.push_front(x);
            } else if (cur_sz < max_sz) {  // add one without popping
                buffer.push_front(x);
                ++cur_sz;
            } else {
                cerr << "Fatal error in LookbackWindow::push()" << endl; exit(1);
            }
        }
        double mean() { return mean(win_sz, front_tm); }
        double mean(size_t lim, size_t ft = 0) {
            if (cur_sz == 0) { cerr << "LookbackWindow.mean(), cur_sz == 0" << endl; exit(1); }
            if (lim + ft > max_sz) {
                cerr << "LookbackWindow::mean(), (lim + ft) > max_sz" << endl;
                exit(1);
            } else if (lim + ft > cur_sz) {
                lim = cur_sz - ft;
            }
            double sum = 0;
            for (size_t i = ft; i < lim + ft; ++i) {
                sum += buffer[i];   // not general for large values, may overflow
            }
            return(sum / lim);
        }
};  // end of class LookbackWindow


//-------------------------------------


bool method_gwa(VcfStripmine::VariantConstPtr_vector& vars);


//-------------------------------------


int 
main(int argc, char* argv[]) {

    VcfStripmine vcfmine;

    processCommandLine(argc, argv, vcfmine);

    // initiate the stripmine, which fetches the first variant from each file
    if (! vcfmine.initiate()) 
        exitusage("Could not initiate stripmine");

    //vcfmine.print();

    // check integrity of sample names
    VcfStripmine::VariantCallFileConstPtr_vector vcfs;
    if (! vcfmine.get_files(vcfs))
        exitmessage("must provide VCF files");
    bool is_fatal = false;
    for (size_t i = 0; i < vcfs.size(); ++i) {
        if (DEBUG(1)) 
            cerr << "main(): " << vcfs[i]->filename << " " << vcfs[i]->sampleNames.size() << " samples:";
        for (size_t j = 0; j < vcfs[i]->sampleNames.size(); ++j) {
            if (vcfs[i]->sampleNames.size() != vcfs[0]->sampleNames.size()) {
                cerr << "*** sample numbers do not match for VCFs: " << vcfs[0]->filename << vcfs[i]->filename << endl;
                is_fatal = true;
            }
            if (vcfs[i]->sampleNames[j] != vcfs[0]->sampleNames[j]) {
                cerr << "*** sample names do not match for VCFs: " << vcfs[0]->filename << vcfs[i]->filename << endl;
                is_fatal = true;
            }
            if (DEBUG(1)) 
                cerr << " " << vcfs[i]->sampleNames[j];
        }
        if (DEBUG(1)) 
            cerr << endl;
    }

    if (is_fatal) 
        exitmessage("sample inconsistencies");

    if (DEBUG(1)) 
        cerr << "main(): All sample names match" << endl;

    VcfStripmine::VariantConstPtr_vector vars;
    while (vcfmine.get(vars)) {

        // for each line of variants in each VCF

        // Variant:: will check for integrity of sample numbers within each VCF
        
        if (samla_method == "gwa") {
            if (! method_gwa(vars))
                exitmessage("failure during method_gwa()");
        } else {
            exitmessage("unknown method '", samla_method, "'");
        }

    }

    return(EXIT_SUCCESS);
}

// small class looking back to previous values, used for calculating contextual quality
static LookbackWindow qwin_Gen(max_gwa_window_size, opt_gwa_window_size, 1);
static LookbackWindow qwin_Wga(max_gwa_window_size, opt_gwa_window_size, 1);
static LookbackWindow qwin_All(max_gwa_window_size, opt_gwa_window_size, 1);


Variant method_gwa_case1(Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case2(Variant& v_Gen, Variant& v_Wga, Variant& v_All, const char* case_name = "2");
Variant method_gwa_case3(Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case4(Variant& v_Gen, Variant& v_Wga, Variant& v_All, const char* case_name = "4");
Variant method_gwa_case5(Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case6(Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case7(Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case8(Variant& v_Gen, Variant& v_Wga, Variant& v_All);


// Implement gwa method, see comments below for details
bool method_gwa(VcfStripmine::VariantConstPtr_vector& vars) {

    bool skip_case = false; // for now

    Variant v_Gen;
    Variant v_Wga;
    Variant v_All;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i]->vcf->filename == vcf_genomic) 
            v_Gen = *vars[i];
        else if (vars[i]->vcf->filename == vcf_wga) 
            v_Wga = *vars[i];
        else if (vars[i]->vcf->filename == vcf_all) 
            v_All = *vars[i];
        else 
            exitmessage("unknown VCF file ", vars[i]->vcf->filename);
    }

    if (DEBUG(2)) {
        if (skip_case && v_Gen.filter == "." && v_Wga.filter == "." && v_All.filter == ".") {  // issue no variant
            return(true);
        }
        if (DEBUG(3)) {
            cerr << "v_Gen: " << v_Gen.vcf->filename << endl;
            cerr << "v_Wga: " << v_Wga.vcf->filename << endl;
            cerr << "v_All: " << v_All.vcf->filename << endl;
        }
    }

    qwin_Gen.push(v_Gen.quality);
    qwin_Wga.push(v_Wga.quality);
    qwin_All.push(v_All.quality);

    /* The possible values for col7 FILTER in these VCFs are:

       .                             matches reference
       LowQual                       low quality but variant seems present
       VQSRTrancheSNP99.00to99.90    low quality just on the edge
       VQSRTrancheSNP99.90to100.00   low quality even juster on the edge
       PASS                          ta-da

       For Genomic, WGA and All libs, if

  ***1  PASS          PASS          x Pass - Combine SNP qualities, get genotypes DP and PL from AllLibs call

     2a PASS          VSQR/LowQual  x Pass - Combine SNP qualities, get genotypes DP and PL from AllLibs call
     2b VSQR/LowQual  PASS          x Pass - Combine SNP qualities, get genotypes DP and PL from AllLibs call
     3  VSQR/LowQual  VSQR/LowQual  x Combine SNP qualities and see if pass certain threshold 
                                      (30 was used in SNPcalling), get genotypes DP and PL from AllLibs call
     4a VSQR/LowQual  No-call       x If no SNP in AllLibs exclude, if SNP in AllLibs also exclude (only 
                                      around 2,000 SNPs are like this)
     4b No-call       VSQR/LowQual  x If no SNP in AllLibs exclude, if SNP in AllLibs also exclude (only
                                      around 2,000 SNPs are like this)
     5  PASS          No-call       x Pass - take  SNP qualities and genotypes from file that it PASSes  (maybe 
                                      there was just no coverage in other - and if it is good enough to Pass in 
                                      one it should be OK)
     6  No-call       PASS          x Pass - take  SNP qualities and genotypes from file that it PASSes  (maybe 
                                      there was just no coverage in other - and if it is good enough to Pass in 
                                      one it should be OK)
  ***7  No-call       No-call       x No call
    */
    
    if (DEBUG(2)) cerr << "**G " << v_Gen << endl << "**W " << v_Wga << endl << "**A " << v_All << endl;

    if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))  // CASE 1
        && 
        (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        Variant v_ans(method_gwa_case1(v_Gen, v_Wga, v_All));
        cout << v_ans << endl;

    } else if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))  // CASE 2a
               && 
               (v_Wga.filter.substr(0, 4) == "VQSR" || v_Wga.filter == "LowQual")) {

        Variant v_ans(method_gwa_case2(v_Gen, v_Wga, v_All, "2a"));
        cout << v_ans << endl;

    } else if ((v_Gen.filter.substr(0, 4) == "VQSR" || v_Gen.filter == "LowQual")  // CASE 2b
               && 
               (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        Variant v_ans(method_gwa_case2(v_Gen, v_Wga, v_All, "2b"));
        cout << v_ans << endl;

    } else if ((v_Gen.filter.substr(0, 4) == "VQSR" || v_Gen.filter == "LowQual")   // CASE 3
               && 
               (v_Wga.filter.substr(0, 4) == "VQSR" || v_Wga.filter == "LowQual")) {

        Variant v_ans(method_gwa_case3(v_Gen, v_Wga, v_All));
        cout << v_ans << endl;

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                     // CASE 4a
               && 
               (v_Wga.filter.substr(0, 4) == "VQSR" || v_Wga.filter == "LowQual")) {

        Variant v_ans(method_gwa_case4(v_Gen, v_Wga, v_All, "4a"));
        cout << v_ans << endl;

    } else if ((v_Gen.filter.substr(0, 4) == "VQSR" || v_Gen.filter == "LowQual")   // CASE 4b
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")) {

        Variant v_ans(method_gwa_case4(v_Gen, v_Wga, v_All, "4b"));
        cout << v_ans << endl;

    } else if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))  // CASE 5
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")) {

        Variant v_ans(method_gwa_case5(v_Gen, v_Wga, v_All));
        cout << v_ans << endl;

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                    // CASE 6
               && 
               (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        Variant v_ans(method_gwa_case6(v_Gen, v_Wga, v_All));
        cout << v_ans << endl;

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                    // CASE 7
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")
               && 
               !skip_case) {

        Variant v_ans(method_gwa_case7(v_Gen, v_Wga, v_All));
        cout << v_ans << endl;

    } else {
        cerr << "=================  Unhandled case in method_gwa()" << endl;
        cerr << "**G " << v_Gen << endl << "**W " << v_Wga << endl << "**A " << v_All << endl;
        cerr << "---" << endl;
        exit(1);
    }
    if (DEBUG(2)) cout << "---" << endl;

    return(true);
}


string
generate_context_qual(const Variant& v_Gen, const Variant& v_Wga, const Variant& v_All) {
    stringstream ss;
    ss << "G:" << showpos << (v_Gen.quality - qwin_Gen.mean()) << ",";
    ss << "W:" << showpos << (v_Wga.quality - qwin_Wga.mean()) << ",";
    ss << "A:" << showpos << (v_All.quality - qwin_All.mean());
    return(ss.str());
}


Variant
method_gwa_case1(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // 1  Unambiguous SNP call, combine strength from v_Gen and v_Wga.  Start with
    // v_All, genotypes and info stats are good, but adjust variant quality to be the
    // sum of v_Gen and v_Wga
    if (DEBUG(2)) cout << "*** case 1 method_gwa_case1()" << endl;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        cerr << "**** strange case 1, v_Gen and v_Wga are PASS but v_All is not" << endl;
        exit(1);
    }
    Variant v_ANS(v_All);
    v_ANS.filter = "SamlaGWA1";
    v_ANS.quality = v_Gen.quality + v_Wga.quality;
    v_ANS.info["Samla"].push_back("GWA1");
    v_ANS.info["Samla"].push_back("Qual_G+W");
    v_ANS.info["Samla"].push_back("Snp_A");
    return(v_ANS);
}


Variant
method_gwa_case2(Variant& v_Gen, Variant& v_Wga, Variant& v_All, const char* case_name) {

    // 2  Unambiguous SNP call for one, filtered for the other. Combine strength from 
    // both, using straight quality for the variant and contextual quality for
    // the filtered one.
    //
    // A corner case could be that v_All doesn't contain a variant for us to
    // pull info from.  When this happens, pull the output from the one that has the
    // SNP call and add the contextual quality from the other...

    if (DEBUG(2)) cout << "*** case 2 method_gwa_case2()" << endl;

    string nm = string("GWA") + case_name;
    Variant v_ANS;
    if (string(case_name) == "2a") { // v_Gen is the variant, v_Wga is filtered
        if (v_All.info["VariantType"][0] == "NO_VARIATION") {
            v_ANS = v_Gen;
            v_ANS.info["Samla"].push_back("Snp_G");
        } else {
            v_ANS = v_All;
            v_ANS.info["Samla"].push_back("Snp_A");
        }
        v_ANS.filter = string("Samla") + nm;
        v_ANS.info["Samla"].push_back(nm);
        v_ANS.quality = v_Gen.quality + abs(v_Wga.quality - qwin_Wga.mean());
        v_ANS.info["Samla"].push_back("Qual_G+WContext");
    } else if (string(case_name) == "2b") {  // v_Wga is the variant, v_Gen is filtered
        if (v_All.info["VariantType"][0] == "NO_VARIATION") {
            v_ANS = v_Wga;
            v_ANS.info["Samla"].push_back("Snp_W");
        } else {
            v_ANS = v_All;
            v_ANS.info["Samla"].push_back("Snp_A");
        }
        v_ANS.filter = string("Samla") + nm;
        v_ANS.info["Samla"].push_back(nm);
        v_ANS.quality = abs(v_Gen.quality - qwin_Gen.mean()) + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_W+GContext");
    } else {
        cerr << "method_gwa_case2(): Unhandled case_name " << case_name << endl;
    }
    v_ANS.info["SamlaContextQual"].push_back(generate_context_qual(v_Gen, v_Wga, v_All));
    if (DEBUG(2) && v_All.info["VariantType"][0] == "NO_VARIATION") {
        cerr << "**** unusual case 2, v_Gen and/or v_Wga are PASS with other filtered but v_All is not a variant" << endl;
        cerr << "**** G " << v_Gen << endl << "**** W " << v_Wga << endl << "**** A " << v_All << endl;
        cerr << "**** S " << v_ANS << endl;
    }
    return(v_ANS);
}


Variant
method_gwa_case3(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // 3  Low-qual genotypes for v_Gen and v_Wga, so do we pass a threshold to emit a genotype?
    if (DEBUG(2)) cout << "*** case 3 method_gwa_case3()" << endl;
    Variant v_ANS;
    if (DEBUG(3)) {
        cerr << "G "; qwin_Gen.dump(10, 1);
        cerr << "W "; qwin_Wga.dump(10, 1);
        cerr << "A "; qwin_All.dump(10, 1);
    }
    double qdelta_Gen = (v_Gen.quality - qwin_Gen.mean());
    double qdelta_Wga = (v_Wga.quality - qwin_Wga.mean());
    double qdelta_All = (v_All.quality - qwin_All.mean());
    double qdelta_Contextual = abs(qdelta_Gen) + abs(qdelta_Wga);
    bool   passedContextualQuality = (qdelta_Contextual >= opt_gwa_quality);
    if (DEBUG(1)) {
        fprintf(stderr, "G quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Gen.quality, qwin_Gen.mean(), qdelta_Gen);
        fprintf(stderr, "W quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Wga.quality, qwin_Wga.mean(), qdelta_Wga);
        fprintf(stderr, "A quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_All.quality, qwin_All.mean(), qdelta_All);
        fprintf(stderr, "                  Contextual abs(G)+abs(W) quality delta: %9.4f ", qdelta_Contextual);
        cerr << (passedContextualQuality ? "PASSED" : "FAILED") << " gwa-quality=" << opt_gwa_quality << endl;
    }
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        // v_All does not have a variant, is this good or bad?
        v_ANS = v_All;
        v_ANS.filter = "SamlaGWA3";
        v_ANS.info["Samla"].push_back("GWA3");
        v_ANS.quality = v_Gen.quality + qdelta_Contextual;
        v_ANS.info["Samla"].push_back("Qual_G+Context");
        v_ANS.info["Samla"].push_back("NoSnp_A");
        // In all cases our quality is v_Gen + qdelta_Contextual
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We think there should be a variant here, because of the drop in
            // contextual quality for both Gen and Wga that passes the quality
            // filter, but we cannot see it because GATK doesn't emit complete
            // information.
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("INFER_UNCALLED_VARIANT");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
        } else {  // mixed bag, either contextual is too low or inconsistent
            v_ANS.filter = "SamlaGWA3";
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
        }
    } else {  
        // v_All has a variant, but we don't yet know if we accept it
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We accept the v_All variant
            v_ANS = v_All;
            v_ANS.filter = "SamlaGWA3";
            v_ANS.info["Samla"].push_back("GWA3");
            v_ANS.quality = v_Gen.quality + qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_G+Context");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("AcceptSnp_A");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            // We can't accept the v_All Variant, build on the v_Gen entry
            v_ANS = v_Gen;
            v_ANS.filter = "SamlaGWA3";
            v_ANS.info["Samla"].push_back("GWA3");
            v_ANS.quality = v_Gen.quality + qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_G+Context");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        } else {
            v_ANS = v_Gen;
            v_ANS.filter = "SamlaGWA3";
            v_ANS.info["Samla"].push_back("GWA3");
            v_ANS.quality = v_Gen.quality + qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_G+Context");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        }
    }
    v_ANS.info["SamlaContextQual"].push_back(generate_context_qual(v_Gen, v_Wga, v_All));
    return(v_ANS);
}


Variant
method_gwa_case4(Variant& v_Gen, Variant& v_Wga, Variant& v_All, const char* case_name) {
    // 4  No SNP in v_Gen and VQSR/LowQual in WGA, emit no SNP (despite what v_All might show)
    // Favour v_Gen if v_All disagrees
    if (DEBUG(2)) cout << "*** case 4 method_gwa_case4()" << endl;
    Variant v_ANS;
    string nm = string("GWA") + case_name;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.filter = string("Samla") + nm;
        v_ANS.info["Samla"].push_back(nm);
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp_A");
    } else {
        v_ANS = v_Gen;
        v_ANS.filter = string("Samla") + nm;
        v_ANS.info["Samla"].push_back(nm);
        v_ANS.info["Samla"].push_back("Qual_G");
        v_ANS.info["Samla"].push_back("Snp_G");
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
    }
    return(v_ANS);
}


Variant
method_gwa_case5(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // 5  SNP in v_Gen and no SNP in v_Wga, accept the v_Gen SNP
    if (DEBUG(2)) cout << "*** case 5 method_gwa_case5()" << endl;
    Variant v_ANS(v_Gen);
    v_ANS.filter = "SamlaGWA5";
    v_ANS.info["Samla"].push_back("Qual_G");
    v_ANS.info["Samla"].push_back("Snp_G");
    return(v_ANS);
}


Variant
method_gwa_case6(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // 6  no SNP in v_Gen and SNP in v_Wga, accept the v_Wga SNP
    if (DEBUG(2)) cout << "*** case 6 method_gwa_case6()" << endl;
    Variant v_ANS(v_Wga);
    v_ANS.filter = "SamlaGWA6";
    v_ANS.info["Samla"].push_back("Qual_W");
    v_ANS.info["Samla"].push_back("Snp_W");
    return(v_ANS);
}


Variant
method_gwa_case7(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // 7  no SNP in v_Gen and no SNP in v_Wga, no SNP
    if (DEBUG(2)) cout << "*** case 7 method_gwa_case7()" << endl;
    if (v_All.info["VariantType"][0] != "NO_VARIATION") {
        cerr << "**** strange case 7, v_Gen and v_Wga are '.' but v_All is *not* NO_VARIATION" << endl;
        exit(1);
    }
    Variant v_ANS(v_All);
    v_ANS.filter = "SamlaGWA7";
    v_ANS.quality = v_Gen.quality + v_Wga.quality;
    v_ANS.info["Samla"].push_back("Qual_G+W");
    v_ANS.info["Samla"].push_back("NoSnp");
    return(v_ANS);
}

