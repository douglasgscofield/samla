// samla.cpp  (c) Douglas G. Scofield, Dept. Plant Physiology, Umeå University
//
// samla [Swedish]: (v) to collect
// Collect VCF files containing variants from multiple sources and interpret
// the results
//
// Uses the vcflib (https://github.com/ekg/vcflib) library for reading VCF
// files.  I am using my own fork (https://github.com/douglasgscofield/vcflib)
// that is changed to:
//   - use const a bit more often especially on const accessors
//   - add Vcf.filename field
//   - TODO rewrite operator<< to not modify fields of Variant
//   - TODO add method Variant.fixup to make modifications to Variant that
//     are currently done by operator<<
//   - TODO handle missing data on output like the VCF presented it on input
//     e.g. if VCF had genotype ./. output this instead of .
//

#define NAME "samla"
#define VERSION "0.0.2"

// CHANGELOG
//
// 0.0.2 : Implement method 'gwa'
//
// 0.0.1 : Implement and lightly test VcfStripmine class

// TODO
// --- Deal with indels in some way
// --- Generalise a method for determining whether variant is no-variant,
//     to bypass gwa reliance on GATK's VariantType=NO_VARIATION.  Perhaps
//     derive a class from Variant adding a method/field to do each variant?
// --- Once the above is done, generalise gwa to not rely on GATK's 
//     VariantType=NO_VARIATION
// --- For gwa method, turn lookback window into surrounding window.  This
//     is a bit tricky as it will involve caching
// --- Return Variants from VcfStripmine in order of VCFs
//

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>

// vcflib, for this to work, requires -I/path/to/vcflib during compilation
//
#include "src/Variant.h"

// very nice argument handling, from https://code.google.com/p/simpleopt/
// TODO update to version 3.6
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
static double       opt_gwa_vqsr_quality = 30.0;
static double       opt_gwa_lowqual_quality = 30.0;
static double       opt_gwa_lowqual_quality_ref = 20.0;
static double       opt_gwa_mixed_quality = 30.0;
static bool         opt_filter_annotate = true;
static bool         opt_stdio = false;
static size_t       opt_debug = 1;
#define             DEBUG(__level__) (opt_debug >= __level__)
static size_t       debug_progress = 10000;
static size_t       opt_progress = 0; // 1000000;
static size_t       num_variants = 0;
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
                        // delete vcf;  // I think this is deleted by delete var below
                        // delete var;  // These statements caused a crash on uppmax
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
            // check integrity of sample names: must match up among all VCFs
            // TODO: do they have to match up when just one sample per file?
            VariantCallFileConstPtr_vector vcfs;
            get_files(vcfs);
            bool is_fatal = false;
            for (size_t i = 0; i < vcfs.size(); ++i) {
                if (DEBUG(1)) 
                    cerr << ".initiate(): " << vcfs[i]->filename << " " << vcfs[i]->sampleNames.size() << " samples:";
                for (size_t j = 0; j < vcfs[i]->sampleNames.size(); ++j) {
                    if (vcfs[i]->sampleNames.size() != vcfs[0]->sampleNames.size()) {
                        cerr << ".initiate(): sample numbers do not match for VCFs: " << vcfs[0]->filename << vcfs[i]->filename << endl;
                        is_fatal = true;
                    }
                    if (vcfs[i]->sampleNames[j] != vcfs[0]->sampleNames[j]) {
                        cerr << ".initiate(): sample names do not match for VCFs: " << vcfs[0]->filename << vcfs[i]->filename << endl;
                        is_fatal = true;
                    }
                    if (DEBUG(1)) 
                        cerr << " " << vcfs[i]->sampleNames[j];
                }
                if (DEBUG(1)) 
                    cerr << endl;
            }
            if (is_fatal) 
                return false;
            if (DEBUG(1)) 
                cerr << ".initiate(): All sample names match" << endl;
            // ready to return
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
           o_filter_annotate, o_no_filter_annotate,
           o_vcf_genomic, o_vcf_wga, o_vcf_all, 
           o_gwa_window, o_gwa_quality, o_gwa_vqsr_quality, o_gwa_lowqual_quality, o_gwa_lowqual_quality_ref, o_gwa_mixed_quality,
           o_output, o_stdio, o_debug, o_progress, o_help };

    CSimpleOpt::SOption smorgas_options[] = {
        { o_references,  "-r",            SO_REQ_SEP },  // file of reference sequence names
        { o_references,  "--references",  SO_REQ_SEP },
        { o_filter_annotate,    "--filter-annotate",    SO_NONE },
        { o_no_filter_annotate, "--no-filter-annotate", SO_NONE },
        { o_method,      "-m",            SO_REQ_SEP },  // method for variant combining
        { o_method,      "--method",      SO_REQ_SEP },

        { o_vcf_genomic, "--vcf-genomic", SO_REQ_SEP },
        { o_vcf_genomic, "-g",            SO_REQ_SEP },
        { o_vcf_wga,     "--vcf-wga",     SO_REQ_SEP },
        { o_vcf_wga,     "-w",            SO_REQ_SEP },
        { o_vcf_all,     "--vcf-all",     SO_REQ_SEP },
        { o_vcf_all,     "-a",            SO_REQ_SEP },
        { o_gwa_window,              "--gwa-window",              SO_REQ_SEP },
        { o_gwa_quality,             "--gwa-quality",             SO_REQ_SEP },
        { o_gwa_vqsr_quality,        "--gwa-vqsr-quality",        SO_REQ_SEP },
        { o_gwa_lowqual_quality,     "--gwa-lowqual-quality",     SO_REQ_SEP },
        { o_gwa_lowqual_quality_ref, "--gwa-lowqual-quality-ref", SO_REQ_SEP },
        { o_gwa_mixed_quality,       "--gwa-mixed-quality",       SO_REQ_SEP },

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
            case o_filter_annotate:    
                opt_filter_annotate = true; break;
            case o_no_filter_annotate:    
                opt_filter_annotate = false; break;
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
                opt_gwa_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_quality; 
                // set all quality threshholds to this value
                opt_gwa_vqsr_quality = opt_gwa_lowqual_quality = opt_gwa_lowqual_quality_ref = opt_gwa_mixed_quality = opt_gwa_quality;
                break;
            case o_gwa_vqsr_quality:   
                opt_gwa_vqsr_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_vqsr_quality; break;
            case o_gwa_lowqual_quality:   
                opt_gwa_lowqual_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_lowqual_quality; break;
            case o_gwa_lowqual_quality_ref:   
                opt_gwa_lowqual_quality_ref = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_lowqual_quality_ref; break;
            case o_gwa_mixed_quality:   
                opt_gwa_mixed_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_mixed_quality; break;
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
    cerr << NAME << " " << VERSION << endl;
    cerr << endl;
    cerr << "Usage:   " << NAME << " [options] -r refnames.txt <in1.vcf> [ in2.vcf ... ]" << endl;
    cerr << endl;
    cerr << "Collect VCF results and produce a consensus list of variants." << endl;
    cerr << endl;
    cerr << "NOTE: " << NAME << " " << VERSION << " is under active development." << endl;
    cerr << endl;
    cerr << "     --references FILE            file containing reference names in order [REQUIRED]" << endl;
    cerr << "                                  should be in the same order as the VCF files" << endl;
    cerr << "     --output FILE                output file name [default is stdout]" << endl;
    cerr << "     --debug INT                  debug info level INT [" << opt_debug << "]" << endl;
    cerr << "     --progress INT               print variants processed mod INT [" << opt_progress << "]" << endl;
    cerr << endl;
    cerr << "     --filter-annotate            annotate FILTER field with additional method-specific information" << endl;
    cerr << "     --no-filter-annotate         do not annotate FILTER field, use only PASS and FAIL" << endl;
    cerr << endl;
    cerr << "     --method METHOD              use combining method 'METHOD', only 'gwa' is implemented" << endl;
    cerr << endl;
    cerr << "'gwa' method options:" <<endl;
    cerr << endl;
    cerr << "     --gwa-window INT                 lookback window size for mean quality, max " << max_gwa_window_size << "[" << opt_gwa_window_size << "]" << endl;
    cerr << "     --gwa-quality FLOAT              minimum quality when combining both VQSR and LowQual variants [" << opt_gwa_quality << "]" << endl;
    cerr << "                                      Specifying this option will set all the quality values below to the given value." << endl;
    cerr << "     --gwa-vqsr-quality FLOAT         minimum quality when combining VQSR variants [" << opt_gwa_vqsr_quality << "]" << endl;
    cerr << "     --gwa-lowqual-quality FLOAT      minimum quality when combining LowQual variants and call does not match reference [" << opt_gwa_lowqual_quality << "]" << endl;
    cerr << "     --gwa-lowqual-quality-ref FLOAT  minimum quality to meet when combining LowQual variants and call matches reference [" << opt_gwa_lowqual_quality_ref << "]" << endl;
    cerr << "     --gwa-mixed-quality FLOAT        minimum quality when combining VQSR with LowQual variants [" << opt_gwa_mixed_quality << "]" << endl;
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

    VcfStripmine::VariantConstPtr_vector vars;

    while (vcfmine.get(vars)) {

        // for each common set of positions containing variants in one or more VCF

        if (samla_method == "gwa") {

            if (! method_gwa(vars))
                exitmessage("failure during method_gwa()");

        } else {

            exitmessage("unknown method '", samla_method, "'");

        }

        ++num_variants;
        if (opt_progress && (num_variants % opt_progress == 0)) {
            cerr << NAME "::main(): processed " << num_variants << " positions, last was " 
                << vars[0]->sequenceName << ":" << vars[0]->position << endl;
        }

    }

    return(EXIT_SUCCESS);
}

// small class looking back to previous values, used for calculating contextual quality
static LookbackWindow qualwindow_Gen(max_gwa_window_size, opt_gwa_window_size, 1);
static LookbackWindow qualwindow_Wga(max_gwa_window_size, opt_gwa_window_size, 1);
static LookbackWindow qualwindow_All(max_gwa_window_size, opt_gwa_window_size, 1);

Variant method_gwa_case1   (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case2_G (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case2_W (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case3   (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case4   (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case5   (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case6_G (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case6_W (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case7_G (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case7_W (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case8_G (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case8_W (Variant& v_Gen, Variant& v_Wga, Variant& v_All);
Variant method_gwa_case9   (Variant& v_Gen, Variant& v_Wga, Variant& v_All);

string
generate_gwa_qual_string(const Variant& v_Gen, const Variant& v_Wga, const Variant& v_All) {
    stringstream ss;
    ss << "G:" << v_Gen.quality << "/" << showpos << (v_Gen.quality - qualwindow_Gen.mean()) << noshowpos << ";";
    ss << "W:" << v_Wga.quality << "/" << showpos << (v_Wga.quality - qualwindow_Wga.mean()) << noshowpos << ";";
    ss << "A:" << v_All.quality << "/" << showpos << (v_All.quality - qualwindow_All.mean()) << noshowpos;
    return(ss.str());
}


void
set_filter(Variant& v, const string& filter) {
    v.filter = filter;
    v.info["Samla"].push_back(filter);
}

void
annotate_filter(Variant& v, const string& filter) {
    if (opt_filter_annotate == true) {
        v.addFilter(string("Samla") + filter);
    }
    v.info["Samla"].push_back(filter);
}

// Implement gwa method, see comments below for details
bool method_gwa(VcfStripmine::VariantConstPtr_vector& vars) {

    bool skip_case = false; // set to true to skip the "easy" cases 9s

    // Make mutable copies of G, W and A variants.  Their order seems to be mixed
    // up in vars[] more than I thought, so be dynamic about determining which is
    // which.
    Variant v_Gen;
    Variant v_Wga;
    Variant v_All;
    for (size_t i = 0; i < vars.size(); ++i) {
        if (vars[i]->vcf->filename == vcf_genomic) {
            v_Gen = *vars[i];
        } else if (vars[i]->vcf->filename == vcf_wga) {
            v_Wga = *vars[i];
        } else if (vars[i]->vcf->filename == vcf_all) {
            v_All = *vars[i];
        } else {
            exitmessage("unknown VCF file ", vars[i]->vcf->filename);
        }
    }
    // Check various features of variants.
    // Note that var.info.size() == 0 will also catch uninitialised vars, though not explicitly.
    stringstream ss;
    if (v_Gen.info.size() == 0)
        ss << "method_gwa(): Genomic variant at " << v_Gen.sequenceName << ":" << v_Gen.position 
            << " has an empty INFO field (col 8), method 'gwa' requires this" << endl;
    if (v_Wga.info.size() == 0)
        ss << "method_gwa(): WGA variant at " << v_Wga.sequenceName << ":" << v_Wga.position 
            << " has an empty INFO field (col 8), method 'gwa' requires this" << endl;
    if (v_All.info.size() == 0)
        ss << "method_gwa(): All variant at " << v_All.sequenceName << ":" << v_All.position 
            << " has an empty INFO field (col 8), method 'gwa' requires this" << endl;
    if (ss.str().length())
        exitmessage(ss.str());

    if (DEBUG(2)) {
        if (skip_case && v_Gen.filter == "." && v_Gen.alleles[1] == "."
                      && v_Wga.filter == "." && v_Wga.alleles[1] == "."
                      && v_All.filter == "." && v_All.alleles[1] == ".") {  // no variants to consider
            return(true);
        }
        if (DEBUG(3)) {
            cerr << "v_Gen: " << v_Gen.vcf->filename << endl;
            cerr << "v_Wga: " << v_Wga.vcf->filename << endl;
            cerr << "v_All: " << v_All.vcf->filename << endl;
        }
    }

    qualwindow_Gen.push(v_Gen.quality);
    qualwindow_Wga.push(v_Wga.quality);
    qualwindow_All.push(v_All.quality);

    /* The possible values for col7 FILTER in these VCFs are:
       .                             matches reference (not always true, this may be set for a called variant too
       LowQual                       low quality but variant seems present
       VQSRTrancheSNP99.00to99.90    low quality just on the edge
       VQSRTrancheSNP99.90to100.00   low quality even juster on the edge
       PASS                          ta-da
    */
    
    if (DEBUG(2)) cerr << "**G " << v_Gen << endl << "**W " << v_Wga << endl << "**A " << v_All << endl;

    Variant v_ANS;

    if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))         // CASE 1
        && 
        (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        v_ANS = method_gwa_case1(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))  // CASE 2_G
               && 
               (v_Wga.filter.substr(0, 4) == "VQSR" || v_Wga.filter == "LowQual")) {

        v_ANS = method_gwa_case2_G(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter.substr(0, 4) == "VQSR" || v_Gen.filter == "LowQual")            // CASE 2_W
               && 
               (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        v_ANS = method_gwa_case2_W(v_Gen, v_Wga, v_All);

    } else if (v_Gen.filter.substr(0, 4) == "VQSR" && v_Wga.filter.substr(0, 4) == "VQSR") { // CASE 3

        v_ANS = method_gwa_case3(v_Gen, v_Wga, v_All);

    } else if (v_Gen.filter == "LowQual" && v_Wga.filter == "LowQual") {                    // CASE 4

        v_ANS = method_gwa_case4(v_Gen, v_Wga, v_All);

// must add 5
    } else if ((v_Gen.filter.substr(0, 4) == "VQSR" && v_Wga.filter == "LowQual")
               ||
               (v_Gen.filter == "LowQual" && v_Wga.filter.substr(0, 4) == "VQSR")) {         // CASE 5

        v_ANS = method_gwa_case5(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                     // CASE 6_G
               && 
               v_Wga.filter.substr(0, 4) == "VQSR") {

        v_ANS = method_gwa_case6_G(v_Gen, v_Wga, v_All);

    } else if (v_Gen.filter.substr(0, 4) == "VQSR"                                  // CASE 6_W
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")) {

        v_ANS = method_gwa_case6_W(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                     // CASE 7_G
               && 
               v_Wga.filter == "LowQual") {

        v_ANS = method_gwa_case7_G(v_Gen, v_Wga, v_All);

    } else if (v_Gen.filter == "LowQual"                                            // CASE 7_W
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")) {

        v_ANS = method_gwa_case7_W(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "PASS" || (v_Gen.filter == "." && v_Gen.alleles[1] != "."))  // CASE 8_G
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")) {

        v_ANS = method_gwa_case8_G(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                    // CASE 8_W
               && 
               (v_Wga.filter == "PASS" || (v_Wga.filter == "." && v_Wga.alleles[1] != "."))) {

        v_ANS = method_gwa_case8_W(v_Gen, v_Wga, v_All);

    } else if ((v_Gen.filter == "." && v_Gen.alleles[1] == ".")                    // CASE 9
               && 
               (v_Wga.filter == "." && v_Wga.alleles[1] == ".")
               && 
               !skip_case) {

        v_ANS = method_gwa_case9(v_Gen, v_Wga, v_All);

    } else {
        cerr << "=================  Unhandled case in method_gwa()" << endl;
        cerr << "**G " << v_Gen << endl << "**W " << v_Wga << endl << "**A " << v_All << endl;
        cerr << "---" << endl;
        exit(1);
    }

    annotate_filter(v_ANS, generate_gwa_qual_string(v_Gen, v_Wga, v_All));

    cout << v_ANS << endl;

    if (DEBUG(2)) cout << "---" << endl;

    return(true);
}

/*********************** Current set of case labels

Case 1: Unambiguous variant for G and W.  What is emitted depends on whether A is
also a variant, and if not which of G and W has the highest quality.  Suffix indicates
which variant was emitted as the variant (_A, _G, _W).

GWA1_A : G, W, A are variants : Emit A as variant. PASS

    The unambiguous variant A is the emitted variant.  Quality is the summed qualities
    of G and W.

GWA1_G : G, W are variants, A is *not* a variant : Emit G as variant. PASS

    Since A was not a variant, we choose one of the others to emit.  Quality G
    >= W so G is the emitted variant.  Quality is quality of G.
    
GWA1_W : G, W are variants, A is *not* a variant : Emit W as variant. PASS

    Since A was not a variant, we choose one of the others to emit.  Quality G
    < W so W is the emitted variant.  Quality is quality of W.
    
****

Case 2: Unambiguous variant for one, filtered variant (both VQSR and LowQual)
for the other.  If A is a variant, it is emitted as the site variant, otherwise
the unambiguous variant is emitted as the site variant, and quality is the sum
of the qualities of the unambiguous variant and the contextual quality of the
other.  Suffix 1 indicates which was the unambiguous variant (_G, _W), suffix 2
indicates which was used for the emitted variant (_A, _G, _W).
TODO: does this quality structure make sense?

GWA_2_G_G : G is variant, W is filtered no-variant, A is no-variant : Emit G as variant. PASS

    The unambiguous variant G is the only variant.  Quality is quality of G.

GWA_2_G_A : G is variant, W is filtered no-variant, A is variant : Emit A as variant. PASS

    Both G and A are variants.  Quality is quality of G plus contextual quality of W.

GWA_2_W_W : G is filtered no-variant, W is variant, A is no-variant : Emit W as variant. PASS

    The unambiguous variant W is the only variant.  Quality is quality of W.

GWA_2_W_A : G is filtered no-variant, W is variant, A is variant : Emit A as variant. PASS

    Both W and A are variants.  Quality is quality of W plus contextual quality of G.

****

Case 3: VQSR no-variants for both v_Gen and v_Wga.  What is emitted depends
on whether the --gwa-vqsr-quality threshold is reached.  Suffix 1 indicates
which source was used for generating the emitted variant or no-variant (_A, _G,
_W), suffix 2 indicates what was the result of the quality threshhold filter.

GWA3_A_CALLED_Pass : G and W VQSR no-variant, A *is* variant : Emit A variant.  PASS

    v_All is variant and the sum of the contextual quality of v_Gen and v_Wga
    is above --gwa-vqsr-quality, so we emit A as the variant.  Quality is the
    sum of contextual qualities for G and W.

GWA3_G_CALLED_FailLowQual : G and W VQSR no-variant, A *is* variant : Emit G as no-variant.  FAIL

    v_All is no-variant and the sum of the contextual quality of v_Gen and
    v_Wga is below --gwa-vqsr-quality, so we do not emit a variant.  The
    quality of G >= W so emit G as the no-variant.  Quality is the sum of
    contextual qualities for G and W.

GWA3_W_CALLED_FailLowQual : G and W VQSR no-variant, A *is* variant : Emit W as no-variant.  FAIL

    v_All is no-variant and the sum of the contextual quality of v_Gen and
    v_Wga is below --gwa-vqsr-quality, so we do not emit a variant.  The
    quality of G < W so emit W as the no-variant.  Quality is the sum of
    contextual qualities for G and W.

GWA3_G_CALLED_FailInconsistent : G, W VQSR no-variant, A *is* variant : Emit G as no-variant.  FAIL

    v_All is variant and the contextual qualities for v_Gen and v_Wga have
    opposite signs, so are not in concert in indicating there is evidence for
    an uncalled variant.  Quality of G >= W so emit G as the no-variant.
    Quality is the sum of contextual qualities for G and W.

GWA3_W_CALLED_FailInconsistent : G, W VQSR no-variant, A *is* variant : Emit W as no-variant.  FAIL

    v_All is variant and the contextual qualities for v_Gen and v_Wga have
    opposite signs, so are not in concert in indicating there is evidence for
    an uncalled variant.  Quality of G < W so emit W as the no-variant.
    Quality is the sum of contextual qualities for G and W.

GWA3_A_UNCALLED_Fail : G, W VQSR no-variant, A no-variant but there should be a variant: Emit A as no-variant.  FAIL

    v_All is no-variant but the sum of the contextual quality of v_Gen and v_Wga
    is above --gwa-vqsr-quality, so we *should* be emitting a variant but we can't
    because v_All doesn't contain any variant information.  Quality is sum of
    contextual qualities for G and W.

GWA3_A_UNCALLED_FailLowQual : G, W VQSR no-variant, A no-variant : Emit no variant.  FAIL

    v_All is no-variant and the sum of the contextual quality of v_Gen and v_Wga
    is below --gwa-vqsr-quality, so we do not emit a variant.  Quality is sum of
    contextual qualities for G and W.

GWA3_A_UNCALLED_FailInconsistent : G, W VQSR no-variant, A no-variant : Emit no variant.  FAIL

    v_All is no-variant and the contextual qualities for v_Gen and v_Wga have
    opposite signs, so are not in concert in indicating there is evidence for
    an uncalled variant.  Quality is sum of contextual qualities for G and W.

****

Case 4: Nearly identical to case 3 but with LowQual no-variants rather than
VQSR no-variants.  LowQual no-variants for both v_Gen and v_Wga.  What is emitted
depends on whether the --gwa-lowqual-quality threshold is reached.  Suffix 1
indicates which source was used for generating the emitted variant or
no-variant (_A, _G, _W), suffix 2 indicates what was the result of the quality
threshhold filter.

Cases below are identical to those for case 3, with GWA4 replacing GWA3, LowQual
replacing VQSR, and --gwa-lowqual-quality replacing --gwa-vqsr-quality.

The exceptions to this is three new cases having to do with testing against
--gwa-lowqual-quality-ref, when we think the position should match the reference.

GWA4_A_CALLED_Pass : G and W LowQual no-variant, A *is* variant : Emit A variant.  PASS

GWA4_G_CALLED_Fail_QualityRef : G, W LowQual no-variant, A *is* variant : Emit G as no-variant.  FAIL

    A is called as a variant but we think this is a mistake, because we fail contextual
    quality but pass G quality + W quality >= gwa-lowqual-quality-ref.  We emit G as the
    no-variant because for quality G >= W, with quality G + W

GWA4_W_CALLED_Fail_QualityRef : G, W LowQual no-variant, A *is* variant : Emit W as no-variant.  FAIL

    A is called as a variant but we think this is a mistake, because we fail contextual
    quality but pass G quality + W quality >= gwa-lowqual-quality-ref.  We emit W as the
    no-variant because for quality G < W, with quality G + W

GWA4_G_CALLED_FailLowQual : G and W LowQual no-variant, A *is* variant : Emit G as no-variant.  FAIL

GWA4_W_CALLED_FailLowQual : G and W LowQual no-variant, A *is* variant : Emit W as no-variant.  FAIL

GWA4_G_CALLED_FailInconsistent : G, W LowQual no-variant, A *is* variant : Emit G as no-variant.  FAIL

GWA4_W_CALLED_FailInconsistent : G, W LowQual no-variant, A *is* variant : Emit W as no-variant.  FAIL

GWA4_A_UNCALLED_Pass_QualityRef : G, W LowQual, A no-variant : Emit A as no-variant.  PASS

    A is not a variant and we fail contextual quality (thus we don't think we
    are missing a variant, as for GWA4_A_UNCALLED_Fail), so we think we actually
    match the reference, as we have enough strength for G quality + W quality 
    >= gwa-lowqual-quality-ref, with quality G + W

GWA4_A_UNCALLED_Fail : G, W LowQual no-variant, A no-variant but there should be a variant: Emit A as no-variant.  FAIL

GWA4_A_UNCALLED_FailLowQual : G, W LowQual no-variant, A no-variant : Emit no variant.  FAIL

GWA4_A_UNCALLED_FailInconsistent : G, W LowQual no-variant, A no-variant : Emit no variant.  FAIL

****

Case 5: Nearly identical to case 3 but with mixed VQSR and LowQual no-variants
rather than just VQSR or LowQual no-variants.  In either case, there are
no-variants for both v_Gen and v_Wga.  What is emitted depends on whether the
--gwa-mixed-quality threshold is reached.  Suffix 1 indicates which source was
used for generating the emitted variant or no-variant (_A, _G, _W), suffix 2
indicates what was the result of the quality threshhold filter.  Because the
quality of VQSR calls tends to be higher, those will be favoured if A is not
a no-variant.

Cases below are identical to those for case 3, with GWA5 replacing GWA3,
'LowQual/VQSR' replacing VQSR, and --gwa-mixed-quality replacing
--gwa-vqsr-quality.

GWA5_A_CALLED_Pass : G and W LowQual/VQSR no-variant, A *is* variant : Emit A variant.  PASS

GWA5_G_CALLED_FailLowQual : G and W LowQual/VQSR no-variant, A *is* variant : Emit G as no-variant.  FAIL

GWA5_W_CALLED_FailLowQual : G and W LowQual/VQSR no-variant, A *is* variant : Emit W as no-variant.  FAIL

GWA5_G_CALLED_FailInconsistent : G, W LowQual/VQSR no-variant, A *is* variant : Emit G as no-variant.  FAIL

GWA5_W_CALLED_FailInconsistent : G, W LowQual/VQSR no-variant, A *is* variant : Emit W as no-variant.  FAIL

GWA5_A_UNCALLED_Fail : G, W LowQual/VQSR no-variant, A no-variant but there should be a variant: Emit A as no-variant.  FAIL

GWA5_A_UNCALLED_FailLowQual : G, W LowQual/VQSR no-variant, A no-variant : Emit no variant.  FAIL

GWA5_A_UNCALLED_FailInconsistent : G, W LowQual/VQSR no-variant, A no-variant : Emit no variant.  FAIL

****

Case 6: No variant in one, VQSR-filtered variant in the other.  Emit no
variant.  Suffix 1 based on which variant was no-variant (_G, _W) and suffix 2
which variant was used for the basis of the emitted no-variant (_A, _G, _W).
 
GWA6_G_A  : G no-variant, W VQSR, A emitted : Emit no variant.  PASS

    v_All is no-variant so use v_All as emitted no-variant.  Quality is
    v_Gen.quality+v_Wga.quality.

GWA6_G_G  : G no-variant, W VQSR, G emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality >= v_Wga quality so use v_Gen as
    emitted no-variant. Quality is v_Gen.quality.

GWA6_G_W  : G no-variant, W VQSR, W emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality < v_Wga quality so use v_Wga as emitted
    no-variant.  Quality is v_Wga.quality.

GWA6_W_A  : G VQSR, W no-variant, A emitted : Emit no variant.  PASS

    v_All is no-variant so use v_All as emitted no-variant.  Quality is
    v_Gen.quality+v_Wga.quality.

GWA6_W_G  : G VQSR, W no-variant, G emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality >= v_Wga quality so use v_Gen as
    emitted no-variant. Quality is v_Gen.quality.

GWA6_W_W  : G VQSR, W no-variant, W emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality < v_Wga quality so use v_Wga as emitted
    no-variant.  Quality is v_Wga.quality.

****

Case 7: No variant in one, LowQual-filtered variant in the other.  Emit no
variant.  As for case 6, suffix 1 based on which variant was no-variant (_G,
_W) and suffix 2 which variant was used for the basis of the emitted no-variant
(_A, _G, _W).
 
GWA7_G_A  : G no-variant, W LowQual, A emitted : Emit no variant.  PASS

    v_All is no-variant so use v_All as emitted no-variant.  Quality is
    v_Gen.quality+v_Wga.quality.

GWA7_G_G  : G no-variant, W LowQual, G emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality >= v_Wga quality so use v_Gen as
    emitted no-variant. Quality is v_Gen.quality.

GWA7_G_W  : G no-variant, W LowQual, W emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality < v_Wga quality so use v_Wga as emitted
    no-variant.  Quality is v_Wga.quality.

GWA7_W_A  : G LowQual, W no-variant, A emitted : Emit no variant.  PASS

    v_All is no-variant so use v_All as emitted no-variant.  Quality is
    v_Gen.quality+v_Wga.quality.

GWA7_W_G  : G LowQual, W no-variant, G emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality >= v_Wga quality so use v_Gen as
    emitted no-variant. Quality is v_Gen.quality.

GWA7_W_W  : G LowQual, W no-variant, W emitted : Emit no variant.  PASS

    v_All *is* a variant, v_Gen quality < v_Wga quality so use v_Wga as emitted
    no-variant.  Quality is v_Wga.quality.

****

Case 8: Variant in one, no-variant in the other.  Emit the one with the variant.
Suffix indicate which was the variant (_G, _W).

GWA8_G  : G variant, W no-variant, G emitted : Emit G variant.  PASS

    v_Gen is emitted variant.  Quality is v_Gen.quality.

GWA8_W  : G no-variant, W LowQual, W emitted : Emit W variant.  PASS

    v_Wga is emitted variant.  Quality is v_Wga.quality.

****

Case 9: No-variant in both.  Emit the highest-quality no-variant.  Suffix
indicates which no-variant was the best (_A, _G, _W).

GWA9_A  : G no-variant, W no-variant, A no-variant, A emitted : Emit A no-variant.  PASS

    v_All is emitted variant.  Quality is v_Gen.quality+v_Wga.quality.

GWA9_G  : G no-variant, W no-variant, A *is* variant, G emitted : Emit G no-variant.  PASS

    v_All *is* a variant, v_Gen quality >= v_Wga quality so v_Gen is emitted
    variant.  Quality is v_Gen.quality.  Probably very rare.

GWA9_W  : G no-variant, W no-variant, A *is* variant, W emitted : Emit W no-variant.  PASS

    v_All *is* a variant, v_Gen quality < v_Wga quality so v_Wga is emitted
    variant.  Quality is v_Wga.quality.  Probably very rare.


 **************************/


// -------- 1  Unambiguous SNP call for both v_Gen and v_Wga


Variant
method_gwa_case1(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // Combine strength from v_Gen and v_Wga.  Start with v_All, genotypes and
    // info stats are good, but adjust variant quality to be the sum of v_Gen
    // and v_Wga
    if (DEBUG(2)) cout << "*** case 1 method_gwa_case1()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] != "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA1_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("Snp_A");
    } else {
        if (v_Gen.quality >= v_Wga.quality) { // favour v_Gen
            v_ANS = v_Gen;  // using quality from v_Gen
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA1_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("Snp_G");
        } else {  // v_Gen.quality < v_Wga.quality, so favour v_Wga
            v_ANS = v_Wga;  // using quality from v_Wga
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA1_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("Snp_W");
        }
        if (DEBUG(2)) {
            cerr << "**1* strange case 1, v_Gen and v_Wga are PASS but v_All is not" << endl;
            cerr << "**1* G " << v_Gen << endl << "**1* W " << v_Wga << endl << "**1* A " << v_All << endl;
        }
    }
    return(v_ANS);
}


// -------- 2_G  Unambiguous SNP call for v_Gen, filtered for v_Wga


Variant
method_gwa_case2_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // Combine strength from both, using straight quality for the variant and
    // contextual quality for the filtered one.
    if (DEBUG(2)) cout << "*** case 2_G method_gwa_case2_G()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_Gen;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA2_G_G");
        v_ANS.info["Samla"].push_back("Snp_G");
        if (DEBUG(2)) {
            cerr << "**2_G* v_Gen is PASS with other filtered but v_All is not a variant" << endl;
            cerr << "**2_G* G " << v_Gen << endl << "**2_G* W " << v_Wga << endl << "**2_G* A " << v_All << endl;
        }
    } else {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA2_G_A");
        v_ANS.info["Samla"].push_back("Snp_A");
    }
    v_ANS.quality = v_Gen.quality + abs(v_Wga.quality - qualwindow_Wga.mean());
    v_ANS.info["Samla"].push_back("Qual_G+WContext");
    return(v_ANS);
}


// -------- 2_W  Unambiguous SNP call for v_Wga, filtered for v_Gen


Variant
method_gwa_case2_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    // Combine strength from both, using straight quality for the variant and
    // contextual quality for the filtered one.
    if (DEBUG(2)) cout << "*** case 2_W method_gwa_case2_W()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_Wga;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA2_W_W");
        v_ANS.info["Samla"].push_back("Snp_W");
        if (DEBUG(2)) {
            cerr << "**2_W** v_Wga is PASS with other filtered but v_All is not a variant" << endl;
            cerr << "**2_W** G " << v_Gen << endl << "**2_W** W " << v_Wga << endl << "**2_W** A " << v_All << endl;
        }
    } else {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA2_W_A");
        v_ANS.info["Samla"].push_back("Snp_A");
    }
    v_ANS.quality = abs(v_Gen.quality - qualwindow_Gen.mean()) + v_Wga.quality;
    v_ANS.info["Samla"].push_back("Qual_W+GContext");
    return(v_ANS);
}


// -------- 3  VQSR genotypes for both v_Gen and v_Wga, so do we pass a threshold to emit a genotype?


Variant
method_gwa_case3(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 3 method_gwa_case3()" << endl;
    Variant v_ANS;
    if (DEBUG(3)) {
        cerr << "G "; qualwindow_Gen.dump(10, 1);
        cerr << "W "; qualwindow_Wga.dump(10, 1);
        cerr << "A "; qualwindow_All.dump(10, 1);
    }
    double qdelta_Gen = (v_Gen.quality - qualwindow_Gen.mean());
    double qdelta_Wga = (v_Wga.quality - qualwindow_Wga.mean());
    double qdelta_All = (v_All.quality - qualwindow_All.mean());
    double qdelta_Contextual = abs(qdelta_Gen) + abs(qdelta_Wga);
    bool   passedContextualQuality = (qdelta_Contextual >= opt_gwa_vqsr_quality);
    if (DEBUG(1)) {
        fprintf(stderr, "3 G quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Gen.quality, qualwindow_Gen.mean(), qdelta_Gen);
        fprintf(stderr, "3 W quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Wga.quality, qualwindow_Wga.mean(), qdelta_Wga);
        fprintf(stderr, "3 A quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_All.quality, qualwindow_All.mean(), qdelta_All);
        fprintf(stderr, "3                   Contextual abs(G)+abs(W) quality delta: %9.4f ", qdelta_Contextual);
        cerr << (passedContextualQuality ? "PASSED" : "FAILED") << " gwa-vqsr-quality=" << opt_gwa_vqsr_quality << endl;
    }
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        // v_All does not have a variant, is this good or bad?
        v_ANS = v_All;
        v_ANS.quality = qdelta_Contextual;
        v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
        v_ANS.info["Samla"].push_back("NoSnp_A");
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We think there should be a variant here, because of the drop in
            // contextual quality for both Gen and Wga that passes the quality
            // filter, but we cannot see it because GATK doesn't emit complete
            // information.
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA3_A_UNCALLED_Fail");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("INFER_UNCALLED_VARIANT");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA3_A_UNCALLED_FailInconsistent");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
        } else {  // mixed bag, either contextual is too low or inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA3_A_UNCALLED_FailLowQual");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
        }
    } else {  
        // v_All has a variant, but we don't yet know if we accept it
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We accept the v_All variant
            v_ANS = v_All;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA3_A_CALLED_Pass");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("AcceptSnp_A");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA3_G_CALLED_FailInconsistent");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA3_W_CALLED_FailInconsistent");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        } else {
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA3_G_CALLED_FailLowQual");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA3_W_CALLED_FailLowQual");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        }
    }
    v_ANS.info["SamlaContextQual"].push_back(generate_gwa_qual_string(v_Gen, v_Wga, v_All));
    return(v_ANS);
}


// -------- 4  LowQual genotypes for both v_Gen and v_Wga, so do we pass a threshold to emit a genotype?
// TODO: combine with case 3 but parameterise it using case name and opt_gwa_*_quality, unless
// we think we will do more specific stuff for cases 3 and/or 4?


Variant
method_gwa_case4(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 4 method_gwa_case4()" << endl;
    Variant v_ANS;
    if (DEBUG(3)) {
        cerr << "G "; qualwindow_Gen.dump(10, 1);
        cerr << "W "; qualwindow_Wga.dump(10, 1);
        cerr << "A "; qualwindow_All.dump(10, 1);
    }
    double qdelta_Gen = (v_Gen.quality - qualwindow_Gen.mean());
    double qdelta_Wga = (v_Wga.quality - qualwindow_Wga.mean());
    double qdelta_All = (v_All.quality - qualwindow_All.mean());
    double qdelta_Contextual = abs(qdelta_Gen) + abs(qdelta_Wga);
    bool   passedContextualQuality = (qdelta_Contextual >= opt_gwa_lowqual_quality);
    bool   passedQualityRef = (v_Gen.quality + v_Wga.quality >= opt_gwa_lowqual_quality_ref);
    if (DEBUG(1)) {
        fprintf(stderr, "4 G quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Gen.quality, qualwindow_Gen.mean(), qdelta_Gen);
        fprintf(stderr, "4 W quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Wga.quality, qualwindow_Wga.mean(), qdelta_Wga);
        fprintf(stderr, "4 A quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_All.quality, qualwindow_All.mean(), qdelta_All);
        fprintf(stderr, "4                   Contextual abs(G)+abs(W) quality delta: %9.4f ", qdelta_Contextual);
        cerr << (passedContextualQuality ? "PASSED" : "FAILED") << " gwa-lowqual-quality=" << opt_gwa_lowqual_quality << endl;
    }
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        // v_All does not have a variant, is this good or bad?
        v_ANS = v_All;
        v_ANS.info["Samla"].push_back("NoSnp_A");
        if (! passedContextualQuality && passedQualityRef) {
            // We do not think we have a variant here, we think we should match the reference
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA4_A_UNCALLED_Pass_QualityRef");
            v_ANS.quality = v_Gen.quality + v_Wga.quality;
            v_ANS.info["Samla"].push_back("Qual_G+W");
            v_ANS.info["Samla"].push_back("PASS_QualityRef");
        } else if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We think there should be a variant here, because of the drop in
            // contextual quality for both Gen and Wga that passes the quality
            // filter, but we cannot see it because GATK doesn't emit complete
            // information.
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA4_A_UNCALLED_Fail");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("INFER_UNCALLED_VARIANT");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA4_A_UNCALLED_FailInconsistent");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
        } else {  // mixed bag, either contextual is too low or inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA4_A_UNCALLED_FailLowQual");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
        }
    } else {  
        // v_All has a variant, but we don't yet know if we accept it
        if (! passedContextualQuality && passedQualityRef) {
            // We think we should actually match the reference
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_G_CALLED_Fail_QualityRef");
                v_ANS.info["Samla"].push_back("NoSnp_G");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_W_CALLED_Fail_QualityRef");
                v_ANS.info["Samla"].push_back("NoSnp_W");
            }
            v_ANS.info["Samla"].push_back("FAIL_QualityRef");
            v_ANS.quality = v_Gen.quality + v_Wga.quality;
            v_ANS.info["Samla"].push_back("Qual_G+W");
            v_ANS.info["Samla"].push_back("RejectSnp_A");
        } else if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We accept the v_All variant
            v_ANS = v_All;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA4_A_CALLED_Pass");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("AcceptSnp_A");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_G_CALLED_FailInconsistent");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_W_CALLED_FailInconsistent");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        } else {
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_G_CALLED_FailLowQual");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA4_W_CALLED_FailLowQual");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        }
    }
    v_ANS.info["SamlaContextQual"].push_back(generate_gwa_qual_string(v_Gen, v_Wga, v_All));
    return(v_ANS);
}


// -------- 5  VQSR filter for one of v_Gen or v_Wga and LowQual for the other, so do we pass a threshold to emit a genotype?


Variant
method_gwa_case5(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 5 method_gwa_case5()" << endl;
    Variant v_ANS;
    if (DEBUG(3)) {
        cerr << "G "; qualwindow_Gen.dump(10, 1);
        cerr << "W "; qualwindow_Wga.dump(10, 1);
        cerr << "A "; qualwindow_All.dump(10, 1);
    }
    double qdelta_Gen = (v_Gen.quality - qualwindow_Gen.mean());
    double qdelta_Wga = (v_Wga.quality - qualwindow_Wga.mean());
    double qdelta_All = (v_All.quality - qualwindow_All.mean());
    double qdelta_Contextual = abs(qdelta_Gen) + abs(qdelta_Wga);
    bool   passedContextualQuality = (qdelta_Contextual >= opt_gwa_mixed_quality);
    if (DEBUG(1)) {
        fprintf(stderr, "5 G quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Gen.quality, qualwindow_Gen.mean(), qdelta_Gen);
        fprintf(stderr, "5 W quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_Wga.quality, qualwindow_Wga.mean(), qdelta_Wga);
        fprintf(stderr, "5 A quality: %9.4f  lookback: %9.4f  quality delta: %9.4f\n", v_All.quality, qualwindow_All.mean(), qdelta_All);
        fprintf(stderr, "5                   Contextual abs(G)+abs(W) quality delta: %9.4f ", qdelta_Contextual);
        cerr << (passedContextualQuality ? "PASSED" : "FAILED") << " gwa-mixed-quality=" << opt_gwa_mixed_quality << endl;
    }
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        // v_All does not have a variant, is this good or bad?
        v_ANS = v_All;
        v_ANS.quality = qdelta_Contextual;
        v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
        v_ANS.info["Samla"].push_back("NoSnp_A");
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We think there should be a variant here, because of the drop in
            // contextual quality for both Gen and Wga that passes the quality
            // filter, but we cannot see it because GATK doesn't emit complete
            // information.
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA5_A_UNCALLED_Fail");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("INFER_UNCALLED_VARIANT");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA5_A_UNCALLED_FailInconsistent");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
        } else {  // mixed bag, either contextual is too low or inconsistent
            set_filter(v_ANS, "FAIL");
            annotate_filter(v_ANS, "GWA5_A_UNCALLED_FailLowQual");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
        }
    } else {  
        // v_All has a variant, but we don't yet know if we accept it
        if (passedContextualQuality && qdelta_Gen < 0 && qdelta_Wga < 0) {
            // We accept the v_All variant
            v_ANS = v_All;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA5_A_CALLED_Pass");
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("PASS_ContextQual");
            v_ANS.info["Samla"].push_back("AcceptSnp_A");
        } else if (passedContextualQuality) {  // contextual is inconsistent
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA5_G_CALLED_FailInconsistent");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA5_W_CALLED_FailInconsistent");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextInconsistent");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        } else {
            // We can't accept the v_All Variant, build on the entry with higher quality
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = v_Gen;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA5_G_CALLED_FailLowQual");
            } else {
                v_ANS = v_Wga;
                set_filter(v_ANS, "FAIL");
                annotate_filter(v_ANS, "GWA5_W_CALLED_FailLowQual");
            }
            v_ANS.quality = qdelta_Contextual;
            v_ANS.info["Samla"].push_back("Qual_GContext+WContext");
            v_ANS.info["Samla"].push_back("Fail_ContextQual");
            v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        }
    }
    v_ANS.info["SamlaContextQual"].push_back(generate_gwa_qual_string(v_Gen, v_Wga, v_All));
    return(v_ANS);
}


// -------- 6_G  No variant in v_Gen, filtered variant (VQSR) in v_Wga.  Emit no variant.


Variant
method_gwa_case6_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 6_G method_gwa_case6_G()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA6_G_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp_A");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        if (v_Gen.quality >= v_Wga.quality) {
            // we will probably most often be here
            v_ANS = v_Gen;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA6_G_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("NoSnp_G");
        } else {
            v_ANS = v_Wga;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA6_G_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("NoSnp_W");
        }
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        if (DEBUG(2)) {
            cerr << "**6_G** v_Gen is '.' with v_Wga LowQual but v_All *is* a variant" << endl;
            cerr << "**6_G** G " << v_Gen << endl << "**6_G** W " << v_Wga << endl << "**6_G** A " << v_All << endl;
        }
    }
    return(v_ANS);
}


// -------- 6_W  Filtered variant (VQSR) in v_Gen, no variant in v_Wga.  Emit no variant.


Variant
method_gwa_case6_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 6_W method_gwa_case6_W()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA6_W_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp_A");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        if (v_Gen.quality >= v_Wga.quality) {
            v_ANS = v_Gen;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA6_W_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("NoSnp_G");
        } else {
            // we will probably most often be here
            v_ANS = v_Wga;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA6_W_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("NoSnp_W");
        }
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        if (DEBUG(2)) {
            cerr << "**6_W** v_Gen is LowQual with v_Wga '.' but v_All *is* a variant" << endl;
            cerr << "**6_W** G " << v_Gen << endl << "**6_W** W " << v_Wga << endl << "**6_W** A " << v_All << endl;
        }
    }
    return(v_ANS);
}


// -------- 7_G  No variant in v_Gen, filtered variant (LowQual) in v_Wga.  Emit no variant.


Variant
method_gwa_case7_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 7_G method_gwa_case7_G()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA7_G_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp_A");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        if (v_Gen.quality >= v_Wga.quality) {
            // we will probably most often be here
            v_ANS = v_Gen;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA7_G_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("NoSnp_G");
        } else {
            v_ANS = v_Wga;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA7_G_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("NoSnp_W");
        }
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        if (DEBUG(2)) {
            cerr << "**7_G** v_Gen is '.' with v_Wga LowQual but v_All *is* a variant" << endl;
            cerr << "**7_G** G " << v_Gen << endl << "**7_G** W " << v_Wga << endl << "**7_G** A " << v_All << endl;
        }
    }
    return(v_ANS);
}


// -------- 7_W  Filtered variant (LowQual) in v_Gen, no variant in v_Wga.  Emit no variant.


Variant
method_gwa_case7_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 7_W method_gwa_case7_W()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA7_W_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp_A");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        if (v_Gen.quality >= v_Wga.quality) {
            v_ANS = v_Gen;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA7_W_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("NoSnp_G");
        } else {
            // we will probably most often be here
            v_ANS = v_Wga;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA7_W_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("NoSnp_W");
        }
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        if (DEBUG(2)) {
            cerr << "**7_W** v_Gen is LowQual with v_Wga '.' but v_All *is* a variant" << endl;
            cerr << "**7_W** G " << v_Gen << endl << "**7_W** W " << v_Wga << endl << "**7_W** A " << v_All << endl;
        }
    }
    return(v_ANS);
}


// -------- 8_G  Variant in v_Gen, no variant in v_Wga.  Emit v_Gen variant.


Variant
method_gwa_case8_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 8_G method_gwa_case8_G()" << endl;
    Variant v_ANS(v_Gen);
    set_filter(v_ANS, "PASS");
    annotate_filter(v_ANS, "GWA8_G");
    v_ANS.info["Samla"].push_back("Qual_G");
    v_ANS.info["Samla"].push_back("Snp_G");
    return(v_ANS);
}


// -------- 8_W  No-variant in v_Gen, variant in v_Wga.  Emit v_Wga variant.


Variant
method_gwa_case8_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 8_W method_gwa_case8_W()" << endl;
    Variant v_ANS(v_Wga);
    set_filter(v_ANS, "PASS");
    annotate_filter(v_ANS, "GWA8_W");
    v_ANS.info["Samla"].push_back("Qual_W");
    v_ANS.info["Samla"].push_back("Snp_W");
    return(v_ANS);
}


// -------- 9  No-variant in v_Gen, no-variant in v_Wga.  Emit no-variant,
//             either A if it is also no-variant or whichever of the others
//             has the highest quality.


Variant
method_gwa_case9(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 9 method_gwa_case9()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION") {
        v_ANS = v_All;
        set_filter(v_ANS, "PASS");
        annotate_filter(v_ANS, "GWA9_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        v_ANS.info["Samla"].push_back("Qual_G+W");
        v_ANS.info["Samla"].push_back("NoSnp");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        // We will probably very rarely end up here
        if (v_Gen.quality >= v_Wga.quality) {
            v_ANS = v_Gen;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA9_G");
            v_ANS.info["Samla"].push_back("Qual_G");
            v_ANS.info["Samla"].push_back("NoSnp_G");
        } else {
            v_ANS = v_Wga;
            set_filter(v_ANS, "PASS");
            annotate_filter(v_ANS, "GWA9_W");
            v_ANS.info["Samla"].push_back("Qual_W");
            v_ANS.info["Samla"].push_back("NoSnp_W");
        }
        v_ANS.info["Samla"].push_back("IgnoreSnp_A");
        if (DEBUG(2)) {
            cerr << "**9** v_Gen and v_Wga both no-variant, but v_All *is* a variant" << endl;
            cerr << "**9** G " << v_Gen << endl << "**9** W " << v_Wga << endl << "**9** A " << v_All << endl;
        }
    }
    return(v_ANS);
}

