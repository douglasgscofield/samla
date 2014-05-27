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

// TODO Must study consistency of context quality
// TODO Standardise filter and INFO fields
// filters must be separated by semicolons, and values do not contain semicolons, 
// commas, whitespace, '=', with comma separating multiple values for a key
// INFO fields must be key=val[,val] with key-value groups separated by semicolons
// key is "Samla<key>"
// things through set_filter() and annotate_filter() are on filter;filter;filter
// set_filter() adds to .filter and to .info["SamlaFilter"]
// annotate_filter() adds to .info["SamlaFilter"] and to .filter on option
// annotate_case() adds the case to "SamlaCase"
// annotate_first_case() adds the case to the first element of "SamlaCase" and to .filter of opt_full_filter_annotate
// annotate_case_add_filter() adds the case to the first element of "SamlaCase" and calls annotate_filter()
// annotate_case_add_full_filter() adds the case to the first element of "SamlaCase" and calls annotate_filter() if opt_full_filter_annotate
//
// things set through v_ANS.info["..."].push_back("...") as INFO: SamlaFilter=filter,filter,filter
// case results are INFO: SamlaCase=val,val,val
//
// VCF annotations added by Samla:
//
// Header:
//
//     ##SamlaVersion=NAME-VERSION
//     ##SamlaReference=http://github.com/douglasgscofield/samla
//
//     ##INFO=<ID=SamlaFilter,Number=.,Type=String,Description="Samla filtering and case decisions">
//     ##INFO=<ID=SamlaCase,Number=.,Type=String,Description="Complete annotation of Samla case">
//
// FILTER field: 
//     PASS or FAIL
// and if not --no-annotate-filter:
//     Samla case applied
//     Samla quality string: "quality:G:G-quality/G-context-quality-signed:W:W-quality/W-context-quality-signed:A:A-quality/A-context-quality-signed"
//
// INFO field:
//     SamlaFilter=
//         PASS or FAIL
//         Samla case applied
//     SamlaCase=
//         Samla case applied
//         Samla quality calculated
//         Variant status of G, W and A
//         Other status
//     culprit:G:string:W:string:A:string
//
// VCF annotations removed by Samla:
// culprit=
// VQSLOD=
// NEGATIVE_TRAIN_SITE
// POSITIVE_TRAIN_SITE
// what else?
// do we ever need to change VariantType?
//
// TODO (working on this): if we pass, just use annotate_case()
//                         if we fail, use annotate_case_add_filter()

#define NAME "samla"
#include "version.h"

// CHANGELOG
//
// 0.0.4 : Cleaning up of produced VCF
//
// 0.0.3 : Reimplementation of cases 3 4 5
//
// 0.0.2 : First implementation of 'gwa'
//
// 0.0.1 : Implement and lightly test VcfStripmine class

// TODO
// --- Deal with indels in some way
// --- Generalise a method for determining whether variant is no-variant,
//     to bypass gwa reliance on GATK's VariantType=NO_VARIATION.  Perhaps
//     derive a class from Variant adding a method/field to do each variant?
// --- Once the above is done, generalise gwa to not rely on GATK's 
//     VariantType=NO_VARIATION
// --- For gwa method, what do we need to modify in the INFO field?  Perhaps
//     remove VariantType=NO_VARIATION for accepted variants is one...
// --- For gwa method, produce headers on output
// --- For gwa method, add appropriate entries to headers
// --- For gwa method, ensure INFO is VCF 4.1-compliant
// --- For gwa method, turn lookback window into surrounding window.  This
//     is a bit tricky as it will involve caching
// --- Return Variants from VcfStripmine in order of VCFs
//

// #define NDEBUG  // uncomment to remove assert() code
#include <assert.h>
#if _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED || _POSIX_C_SOURCE >= 200112L || __DARWIN_C_LEVEL >= 200112L
// for readlink(), part of using /proc/self/exe on Linux see below
#include <unistd.h>  
#endif

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

#define BUFMAX 2048

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
static double       opt_gwa_quality_ref = 20.0;
static double       opt_gwa_vqsr_quality = 30.0;
static double       opt_gwa_lowqual_quality = 30.0;
static double       opt_gwa_lowqual_quality_ref = 20.0;
static double       opt_gwa_mixed_quality = 30.0;
static double       opt_gwa_mixed_quality_ref = 20.0;
static bool         opt_gwa_force_consistency = true;

#define apply_gwa_quality() \
    { \
        opt_gwa_vqsr_quality = opt_gwa_lowqual_quality = opt_gwa_lowqual_quality_ref = \
            opt_gwa_mixed_quality = opt_gwa_mixed_quality_ref = opt_gwa_quality; \
    }

#define apply_gwa_quality_ref() \
    { \
        opt_gwa_lowqual_quality_ref = opt_gwa_mixed_quality_ref = opt_gwa_quality_ref; \
    }

static bool         opt_gwa_enable_context_quality = false;
static bool         opt_gwa_vqsr_vqsr_normal = false;
static bool         opt_gwa_lowqual_lowqual_normal = true;
static bool         opt_gwa_mixed_normal = false;
static bool         opt_filter_annotate = false;
static bool         opt_full_filter_annotate = false;
static bool         opt_stdio = false;
static size_t       opt_debug = 0;
#define             DEBUG(__level__) (opt_debug >= __level__)
static size_t       debug_progress = 10000;
static size_t       opt_progress = 0; // 1000000;
static size_t       num_variants = 0;
static const char tab = '\t';
static const string sep = "\t";
static string       samla_name;   // filled by whatever argv[0] is
static string       given_command_line;  // filled by what the user provided
static string       full_command_line;  // filled as if the user was explicit about every option


//----- Utility function declarations

static void 
exitmessage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "", const string& msg5 = "");

static void 
exitusage(const string& msg = "", const string& msg2 = "", const string& msg3 = "", const string& msg4 = "", const string& msg5 = "");

template<class T>
ostream& 
operator<<(ostream& out, const vector<T>& v);

void 
ToLower(string& s);


//----- CLASS VcfStripmine

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
                string            vcf_filename;
                bool              is_open;
                bool              is_allocated;

                VcfStrip(void) : is_open(false) { }

                ~VcfStrip(void) {
                    if (is_allocated) {
                        // delete vcf;  // I think this is deleted by delete var below
                        // delete var;  // These statements caused a crash on uppmax
                    }
                }

                bool get_next_variant(void) {
                    // Call vcflib's getNextVariant() within this wrapper.  Currently it:
                    // 1. skips indels
                    // 2. sets is_open = false for the strip if no more variants
                    while (vcf->getNextVariant(*var)) {
                        if (VcfStripmine::is_indel(*var)) {
                            if (DEBUG(0)) 
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
            exit(1);  // so not *technically* a const method
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
            // or: the_mine[filename].var = new Variant; the_mine[filename].var->setVariantCallFile(the_mine[filename].vcf);
            the_mine[filename].is_allocated = true;
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
        VariantCallFile* get_VariantCallFilePtr(const string& filename) {
            if (the_mine.size() == 0)
                error_exit(".get_file: no vcf files add()ed");
            if (! the_mine.count(filename))
                error_exit(".get_file: no vcf files add()ed named ", filename);
            return the_mine[filename].vcf;
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


string
generate_given_command_line_string(int argc, char* argv[]) {
    // generate only from input args
    stringstream ss;
    ss << argv[0];
    for (int i = 1; i < argc; ++i) {
        ss << " " << argv[i];
    }
    return ss.str();
}


string
generate_full_command_line_string() {
    // generate only from argv[0] and option values
    stringstream ss;

#if _BSD_SOURCE || _XOPEN_SOURCE >= 500 || _XOPEN_SOURCE && _XOPEN_SOURCE_EXTENDED || _POSIX_C_SOURCE >= 200112L || __DARWIN_C_LEVEL >= 200112L
    // Do this to resolve samla_name more completely, from
    // http://stackoverflow.com/questions/12254980/how-to-get-current-running-file-name-c
    //
    // TODO: more complete
    // http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
    // Mac OS X: _NSGetExecutablePath() (man 3 dyld)
    // Linux: readlink /proc/self/exe
    // Solaris: getexecname()
    // FreeBSD: sysctl CTL_KERN KERN_PROC KERN_PROC_PATHNAME -1
    // FreeBSD if it has procfs: readlink /proc/curproc/file (FreeBSD doesn't have procfs by default)
    // NetBSD: readlink /proc/curproc/exe
    // DragonFly BSD: readlink /proc/curproc/file
    // Windows: GetModuleFileName() with hModule = NULL
    //
    // maybe use cmake instead?
    // 
    char buf[BUFMAX + 1];
    ssize_t l = readlink("/proc/self/exe", buf, BUFMAX);
    if (l > 0) {
        buf[l] = '\0';
        ss << buf;
    } else {  // maybe there wasn't /proc/self/exe
        ss << samla_name;
    }
#else
    ss << samla_name;
#endif
    ss << " " << "--references"  << " " << references_file;
    ss << " " << (opt_full_filter_annotate ? "--full-filter-annotate" : (opt_filter_annotate ? "--filter-annotate" : "--no-filter-annotate"));
    ss << " " << "--method"      << " " << samla_method;
    ss << " " << "--vcf-genomic" << " " << vcf_genomic;
    ss << " " << "--vcf-wga"     << " " << vcf_wga;
    ss << " " << "--vcf-all"     << " " << vcf_all;
    ss << " " << "--gwa-window"  << " " << opt_gwa_window_size;
    // ss << " " << "--gwa-quality"  << " " << opt_gwa_quality;  // this option simply sets several others
    ss << " " << "--gwa-vqsr-quality"  << " " << opt_gwa_vqsr_quality;
    ss << " " << "--gwa-lowqual-quality"  << " " << opt_gwa_lowqual_quality;
    ss << " " << "--gwa-lowqual-quality-ref"  << " " << opt_gwa_lowqual_quality_ref;
    ss << " " << "--gwa-mixed-quality"  << " " << opt_gwa_mixed_quality;
    ss << " " << "--gwa-mixed-quality-ref"  << " " << opt_gwa_mixed_quality_ref;
    ss << " " << (opt_gwa_force_consistency ? "--gwa-force-consistency" : "--gwa-no-force-consistency");
    ss << " " << (opt_gwa_enable_context_quality ? "--gwa-enable-context-quality" : "--gwa-disable-context-quality");
    ss << " " << (opt_gwa_vqsr_vqsr_normal ? "--gwa-vqsr-vqsr-normal" : "--gwa-vqsr-vqsr-fail");
    ss << " " << (opt_gwa_lowqual_lowqual_normal ? "--gwa-lowqual-lowqual-normal" : "--gwa-lowqual-lowqual-fail");
    ss << " " << (opt_gwa_mixed_normal ? "--gwa-mixed-normal" : "--gwa-mixed-fail");
    ss << " " << "--output" << " " << output_file;
    ss << " " << "--debug" << " " << opt_debug;
    ss << " " << "--progress" << " " << opt_progress;
    return ss.str();
}


void
processCommandLine(int argc, char* argv[], VcfStripmine& stripmine) {

    enum { o_references, o_method, 
           o_filter_annotate, o_no_filter_annotate, o_full_filter_annotate,
           o_vcf_genomic, o_vcf_wga, o_vcf_all, 
           o_gwa_window, o_gwa_quality, o_gwa_quality_ref, o_gwa_vqsr_quality, 
           o_gwa_lowqual_quality, o_gwa_lowqual_quality_ref, 
           o_gwa_mixed_quality, o_gwa_mixed_quality_ref,
           o_gwa_force_consistency, o_gwa_no_force_consistency,
           o_gwa_enable_context_quality, o_gwa_disable_context_quality,
           o_gwa_vqsr_vqsr_fail, o_gwa_vqsr_vqsr_normal,
           o_gwa_lowqual_lowqual_fail, o_gwa_lowqual_lowqual_normal,
           o_gwa_mixed_fail, o_gwa_mixed_normal,
           o_output, o_stdio, o_debug, o_progress, o_help };

    CSimpleOpt::SOption samla_options[] = {
        { o_references,  "-r",            SO_REQ_SEP },  // file of reference sequence names
        { o_references,  "--references",  SO_REQ_SEP },
        { o_no_filter_annotate,   "--no-filter-annotate",   SO_NONE },
        { o_filter_annotate,      "--filter-annotate",      SO_NONE },
        { o_full_filter_annotate, "--full-filter-annotate", SO_NONE },
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
        { o_gwa_quality_ref,         "--gwa-quality-ref",         SO_REQ_SEP },
        { o_gwa_vqsr_quality,        "--gwa-vqsr-quality",        SO_REQ_SEP },
        { o_gwa_lowqual_quality,     "--gwa-lowqual-quality",     SO_REQ_SEP },
        { o_gwa_lowqual_quality_ref, "--gwa-lowqual-quality-ref", SO_REQ_SEP },
        { o_gwa_mixed_quality,       "--gwa-mixed-quality",       SO_REQ_SEP },
        { o_gwa_mixed_quality_ref,   "--gwa-mixed-quality-ref",   SO_REQ_SEP },
        { o_gwa_force_consistency,       "--gwa-force-consistency",       SO_NONE },
        { o_gwa_no_force_consistency,    "--gwa-no-force-consistency",    SO_NONE },
        { o_gwa_enable_context_quality,  "--gwa-enable-context-quality",  SO_NONE },
        { o_gwa_disable_context_quality, "--gwa-disable-context-quality", SO_NONE },
        { o_gwa_vqsr_vqsr_fail,          "--gwa-vqsr-vqsr-fail",          SO_NONE },
        { o_gwa_vqsr_vqsr_normal,        "--gwa-vqsr-vqsr-normal",        SO_NONE },
        { o_gwa_lowqual_lowqual_fail,    "--gwa-lowqual-lowqual-fail",    SO_NONE },
        { o_gwa_lowqual_lowqual_normal,  "--gwa-lowqual-lowqual-normal",  SO_NONE },
        { o_gwa_mixed_fail,              "--gwa-mixed-fail",              SO_NONE },
        { o_gwa_mixed_normal,            "--gwa-mixed-normal",            SO_NONE },
        { o_gwa_mixed_fail,              "--gwa-vqsr-lowqual-fail",       SO_NONE },
        { o_gwa_mixed_normal,            "--gwa-vqsr-lowqual-normal",     SO_NONE },

        { o_output,      "-o",            SO_REQ_SEP },  // output file
        { o_output,      "--output",      SO_REQ_SEP },
        { o_stdio,       "-",             SO_NONE }, 
        { o_debug,       "--debug",       SO_REQ_SEP },
        { o_progress,    "--progress",    SO_REQ_SEP },
        { o_help,        "--help",        SO_NONE },
        { o_help,        "-h",            SO_NONE },
        { o_help,        "-?",            SO_NONE }, 
        SO_END_OF_OPTIONS
    };

    samla_name = argv[0];

    given_command_line = generate_given_command_line_string(argc, argv);

    CSimpleOpt args(argc, argv, samla_options);

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
                opt_filter_annotate = true; opt_full_filter_annotate = false; break;
            case o_full_filter_annotate:    
                opt_filter_annotate = true; opt_full_filter_annotate = true; break;
            case o_no_filter_annotate:    
                opt_filter_annotate = false; opt_full_filter_annotate = false; break;
            case o_method:   
                samla_method = args.OptionArg(); 
                ToLower(samla_method);
                if (samla_method == "gwa-ksp") {
                    samla_method = "gwa";
                    opt_gwa_quality = 30.0;
                    apply_gwa_quality();
                    opt_gwa_quality_ref = 20.0;
                    apply_gwa_quality_ref();
                    opt_gwa_force_consistency = true;
                    opt_gwa_enable_context_quality = false;
                    opt_gwa_vqsr_vqsr_normal = false;
                    opt_gwa_mixed_normal = false;
                }
                break;
            case o_vcf_genomic:   
                vcf_genomic = args.OptionArg(); break;
            case o_vcf_wga:   
                vcf_wga = args.OptionArg(); break;
            case o_vcf_all:   
                vcf_all = args.OptionArg(); break;
            case o_gwa_window:   
                opt_gwa_window_size = args.OptionArg() ? atoi(args.OptionArg()) : opt_gwa_window_size; break;

            case o_gwa_quality:   
                if (args.OptionArg()) {
                    opt_gwa_quality = atof(args.OptionArg());
                    apply_gwa_quality();
                } else {
                    exitmessage("error interpreting option value for --gwa-quality '", args.OptionArg(), "'");
                }
                break;

            case o_gwa_quality_ref:   
                if (args.OptionArg()) {
                    opt_gwa_quality_ref = atof(args.OptionArg());
                    apply_gwa_quality_ref();
                } else {
                    exitmessage("error interpreting option value for --gwa-quality-ref '", args.OptionArg(), "'");
                }
                break;

            case o_gwa_vqsr_quality:   
                opt_gwa_vqsr_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_vqsr_quality; break;
            case o_gwa_lowqual_quality:   
                opt_gwa_lowqual_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_lowqual_quality; break;
            case o_gwa_lowqual_quality_ref:   
                opt_gwa_lowqual_quality_ref = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_lowqual_quality_ref; break;
            case o_gwa_mixed_quality:   
                opt_gwa_mixed_quality = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_mixed_quality; break;
            case o_gwa_mixed_quality_ref:   
                opt_gwa_mixed_quality_ref = args.OptionArg() ? atof(args.OptionArg()) : opt_gwa_mixed_quality_ref; break;

            case o_gwa_force_consistency:
                opt_gwa_force_consistency = true; break;
            case o_gwa_no_force_consistency:
                opt_gwa_force_consistency = false; break;

            case o_gwa_disable_context_quality:
                opt_gwa_enable_context_quality = false; break;
            case o_gwa_enable_context_quality:
                opt_gwa_enable_context_quality = true; break;

            case o_gwa_vqsr_vqsr_fail:
                opt_gwa_vqsr_vqsr_normal = false; break;
            case o_gwa_vqsr_vqsr_normal:
                opt_gwa_vqsr_vqsr_normal = true; break;
            case o_gwa_lowqual_lowqual_fail:
                opt_gwa_lowqual_lowqual_normal = false; break;
            case o_gwa_lowqual_lowqual_normal:
                opt_gwa_lowqual_lowqual_normal = true; break;
            case o_gwa_mixed_fail:
                opt_gwa_mixed_normal = false; break;
            case o_gwa_mixed_normal:
                opt_gwa_mixed_normal = true; break;

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

    // handle methods, note 'gwa-ksp' is handled during options processing

    if (samla_method == "gwa") {
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

    full_command_line = generate_full_command_line_string();

}


//----- Utility function definitions

template<class T> ostream& 
operator<<(ostream& out, const vector<T>& v) {
    // ------ print vector as [i]val[,[i+1]val ...]
    // ------ e.g., [0]firstval,[1]secondval,[2]thirdval
    size_t i = 0;
    for (typename vector<T>::const_iterator p = v.begin(); p != v.end(); ++i) {
        out << "[" << i << "]" << *p;
        if (++p != v.end()) out << ",";
    }
    return(out);
}

void 
ToLower(string& s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

static void
exitmessage(const string& msg, const string& msg2, const string& msg3, const string& msg4, const string& msg5) {
    if (! msg.empty()) cerr << endl << "*** " << msg << msg2 << msg3 << msg4 << msg5 << endl;
    exit(1);
}

static void
exitusage(const string& msg, const string& msg2, const string& msg3, const string& msg4, const string& msg5) {
#define show_set(__f__) if (__f__) cerr << " [set]";
    if (! msg.empty()) cerr << endl << "*** " << msg << msg2 << msg3 << msg4 << msg5 << endl;
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
    cerr << "     --filter-annotate            annotate FILTER field with additional method-specific information"; show_set(opt_filter_annotate); cerr << endl;
    cerr << "     --full-filter-annotate       annotate FILTER field with even more method-specific information"; show_set(opt_full_filter_annotate); cerr << endl;
    cerr << "     --no-filter-annotate         do not annotate FILTER field, use only PASS and filters indicating failure described in the VCF header"; show_set(!opt_filter_annotate); cerr << endl;
    cerr << endl;
    cerr << "     --method METHOD              use combining method 'METHOD', only 'gwa' and 'gwa-ksp' are implemented" << endl;
    cerr << endl;
    cerr << "'gwa' method options:" <<endl;
    cerr << endl;
    cerr << "     --gwa-window INT                 lookback window size for mean quality, max " << max_gwa_window_size << "[" << opt_gwa_window_size << "]" << endl;
    cerr << "     --gwa-quality FLOAT              minimum quality when combining both VQSR and LowQual variants [" << opt_gwa_quality << "]" << endl;
    cerr << "                                      Specifying this option will set all the quality values below to the given value." << endl;
    cerr << "     --gwa-quality-ref FLOAT          minimum quality when combining variants and call matches reference [" << opt_gwa_quality_ref << "]" << endl;
    cerr << "                                      Specifying this option will also set --gwa-lowqual-quality-ref and --gwa-mixed-quality-ref." << endl;
    cerr << "     --gwa-vqsr-quality FLOAT         minimum quality when combining VQSR variants [" << opt_gwa_vqsr_quality << "]" << endl;
    cerr << "     --gwa-lowqual-quality FLOAT      minimum quality when combining LowQual variants and call does not match reference [" << opt_gwa_lowqual_quality << "]" << endl;
    cerr << "     --gwa-lowqual-quality-ref FLOAT  minimum quality to meet when combining LowQual variants and call matches reference [" << opt_gwa_lowqual_quality_ref << "]" << endl;
    cerr << "     --gwa-mixed-quality FLOAT        minimum quality when combining VQSR with LowQual variants [" << opt_gwa_mixed_quality << "]" << endl;
    cerr << "     --gwa-mixed-quality-ref FLOAT    minimum quality to meet when combining VQSR with LowQual variants and call matches reference [" << opt_gwa_mixed_quality_ref << "]" << endl;
    cerr << endl;
    cerr << "     --gwa-force-consistency          require G, W and A to agree on variant/no-variant for potentially ambiguous cases 4 and 5"; show_set(opt_gwa_force_consistency); cerr << endl;
    cerr << "     --gwa-no-force-consistency       do not require G, W and A to agree on variant/no-variant for potentially ambiguous cases 4 and 5"; show_set(! opt_gwa_force_consistency); cerr << endl;
    cerr << endl;
    cerr << "     --gwa-disable-context-quality    disables usage of context quality, qualities instead compared directly"; show_set(! opt_gwa_enable_context_quality); cerr << endl;
    cerr << "     --gwa-enable-context-quality     enables usage of context quality"; show_set(opt_gwa_enable_context_quality); cerr << endl;
    cerr << "                                      NOTE: context quality is very experimental and may be incorrect, use at you rown risk" << endl;
    cerr << endl;
    cerr << "     --gwa-vqsr-vqsr-fail             mark as FAIL all cases having both genomic and WGA VQSR-filtered variants"; show_set(! opt_gwa_vqsr_vqsr_normal); cerr << endl;
    cerr << "     --gwa-vqsr-vqsr-normal           PASS/FAIL cases having both genomic and WGA VQSR-filtered variants depending on culprits"; show_set(opt_gwa_vqsr_vqsr_normal); cerr << endl;
    cerr << "     --gwa-lowqual-lowqual-fail       mark as FAIL all cases having both genomic and WGA LowQual-filtered variants"; show_set(! opt_gwa_lowqual_lowqual_normal); cerr << endl;
    cerr << "     --gwa-lowqual-lowqual-normal     PASS/FAIL cases having both genomic and WGA LowQual-filtered variants depending on culprits"; show_set(opt_gwa_lowqual_lowqual_normal); cerr << endl;
    cerr << "     --gwa-mixed-fail                 mark as FAIL all cases having both a VQSR-filtered and a LowQual-filtered variant"; show_set(! opt_gwa_mixed_normal); cerr << endl;
    cerr << "     --gwa-mixed-normal               PASS/FAIL cases having both a VQSR-filtered and a LowQual-filtered variant"; show_set(opt_gwa_mixed_normal); cerr << endl;
    cerr << endl;
    cerr << "                                      --gwa-vqsr-lowqual-fail and --gwa-vqsr-lowqual-normal are synonyms for the above --gwa-mixed-* options" << endl;
    cerr << endl;
    cerr << endl;
    cerr << "The 'gwa-ksp' method sets the following options. Options appearing after this may make further changes to option values." <<endl;
    cerr << endl;
    cerr << "     --method gwa" << endl;
    cerr << "     --gwa-quality 30" << endl;
    cerr << "     --gwa-quality-ref 20" << endl;
    cerr << "     --gwa-force-consistency" << endl;
    cerr << "     --gwa-disable-context-quality" << endl;
    cerr << "     --gwa-vqsr-vqsr-fail" << endl;
    cerr << "     --gwa-mixed-fail" << endl;
    cerr << endl;
    cerr << endl;
    cerr << "For methods 'gwa' and 'gwa-ksp', all VCF files must be specified using these options:" << endl;
    cerr << endl;
    cerr << "     --vcf-genomic FILE           VCF file containing genomic calls" << endl;
    cerr << "     --vcf-wga FILE               VCF file containing whole-genome-amplified calls" << endl;
    cerr << "     --vcf-all FILE               VCF file containing all (pooled) calls" << endl;
    cerr << endl;
    cerr << endl;
    cerr << "     -h | -? | --help             help" << endl;
    cerr << endl;
    cerr << "Version:     " << NAME << " " << VERSION << endl;
    cerr << endl;
    cerr << "Build flags: " << CXXFLAGS << endl;
    cerr << endl;
    exit(1);
}


//------- CLASS LookbackWindow


class LookbackWindow {  // a sliding window that holds values at the current and previous positions
    private:
        size_t        max_sz;    // maximum number of values to hold
        size_t        win_sz;    // window size of values to consider for e.g. mean (<= max_size)
        size_t        cur_sz;    // current number of values held, may be less than max_sz if just cleared.
        size_t        front_tm;  // values to ignore at the front, e.g. to exclude current position
        deque<double> circular;    // circular buffer holding values
    public:
        LookbackWindow(size_t msz = max_gwa_window_size, size_t wsz = opt_gwa_window_size, size_t ft = 0) 
            : max_sz(msz), win_sz(wsz), cur_sz(0), front_tm(ft) { circular.resize(max_sz); };
        size_t cur_size(void) { return (cur_sz); }
        size_t win_size(void) { return (win_sz); }
        size_t max_size(void) { return (max_sz); }
        size_t front_trim(void) { return (front_tm); }
        void   clear(void) { cur_sz = 0; circular.clear(); circular.resize(max_sz); }
        void   dump() { dump(win_sz, front_tm); }
        void   dump(size_t n, size_t ft) { 
            if ((n + ft) > cur_sz) n = cur_sz - ft;
            cerr << "LookbackWindow(max_sz=" << max_sz << ",win_sz=" << win_sz << ",front_tm=" << front_tm << "): ";
            cerr << "["; for (size_t i = 0; i < front_tm; ++i) { cerr << circular[i] << ", "; } cerr << "trimmed] ";
            for (size_t i = ft; i < (n + ft); ++i) {
                cerr << circular[i]; if (i == win_sz + ft) cerr << "|win_sz"; if (i + 1 < n + ft) cerr << ", ";
            }
            cerr << " + " << (cur_sz - n - ft) << " additional values not within " << n + ft << " of the start" << endl;
        }
        void   push(double x) {  // add a value
            if (cur_sz == max_sz) { // pop one off the end before adding
                circular.pop_back();
                circular.push_front(x);
            } else if (cur_sz < max_sz) {  // add one without popping
                circular.push_front(x);
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
                sum += circular[i];   // overflow highly unlikely
            }
            return(sum / lim);
        }
};  // end of class LookbackWindow


//-------------------------------------


bool method_gwa(Variant& v_ANS, VcfStripmine::VariantConstPtr_vector& vars);


//-------------------------------------


enum Header_t { Header_fileformat, Header_FILTER, Header_FORMAT, Header_INFO, Header_contig, Header_reference, Header_misc };
static const char* const HeaderTag[] = { "fileformat", "FILTER", "FORMAT", "INFO", "contig", "reference", "" };
class HeaderAnnotation {
    private:
        class AnnotationTag {
            public:
                Header_t ht;
                string tag, description;
                bool   is_filter;
                AnnotationTag(Header_t _ht, const string& _t, const bool _if, const string& _d = "") 
                    : ht(_ht), tag(_t), description(_d), is_filter(_if) { }
                string str() const {
                    stringstream ss;
                    ss << "##"; 
                    switch(ht) {
                        case Header_fileformat:
                            ss << HeaderTag[ht] << "=" << tag; break;
                        case Header_FILTER:
                            ss << HeaderTag[ht] << "=<ID=" << tag << ",Description=\"" << description << "\">"; break;
                        case Header_FORMAT:
                            ss << HeaderTag[ht] << "=<ID=" << tag << description << ">"; break;
                        case Header_INFO:
                            ss << HeaderTag[ht] << "=<ID=" << tag << description << ">"; break;
                        case Header_contig:
                            ss << HeaderTag[ht] << "=<ID=" << tag << ",length=" << description << ">"; break;
                        case Header_misc:
                            ss << tag << description; break;
                        default:
                            exitmessage("unknown case in HeaderAnnotation::AnnotationTag::str()"); break;
                    }
                    return ss.str();
                }
        };
        vector<AnnotationTag> tags;
        typedef vector<AnnotationTag>::const_iterator vATCI;
    public:
        void add(const Header_t ht, const string& t, const bool is_filter, const string& description = "") {
            if (exists(t))
                exitmessage("FILTERAnnotation::add(): filter ", t, " already exists");
            tags.push_back(AnnotationTag(ht, t, is_filter, description));
        }
        bool exists(const string& tag_for_search) {
            // special cases for opt_full_filter_annotate
            if (opt_full_filter_annotate) {
                if (tag_for_search == "quality:") 
                    return false;
                else if (tag_for_search.substr(0, 8) == "quality:") 
                    return true;
                if (tag_for_search == "culprits:") 
                    return false;
                else if (tag_for_search.substr(0, 8) == "culprits:") 
                    return true;
            }
            // look amongst registered Headers
            for (vATCI t = tags.begin(); t != tags.end(); ++t) 
                if (t->tag == tag_for_search) 
                    return true;
            return false;
        }
        void fill_vcf_header(VariantCallFile * vcf, bool print_full = false) const {
            for (vATCI t = tags.begin(); t != tags.end(); ++t)
                if (print_full || t->is_filter) 
                    vcf->addHeaderLine(t->str());
        }
};


static HeaderAnnotation Headers;


//-------------------------------------


void
prepareVcfHeader(VariantCallFile * header_vcf) {

    header_vcf->removeHeaderLines("##FILTER=<ID=");
    header_vcf->removeHeaderLines("##INFO=<ID=culprit=");

    stringstream line, ss;

    line.str("");
    line << "##source=" << NAME << "-" << VERSION;  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaReference=" << REFERENCE;  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaVersion=\"" << VERSION << "\"";  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##VcflibVersion=\"" << VCFLIB_VERSION << "\"";  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##VcflibReference=\"" << VCFLIB_REFERENCE << "\"";  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaCXXBuildVersion=\"" << CXX_VERSION << "\"";  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaCXXBuildFlags=\"" << CXXFLAGS << "\"";  // provided in version.h by build
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaGivenCommandLine=\"" << given_command_line << "\"";  // assembled here
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##SamlaFullCommandLine=\"" << full_command_line << "\"";  // assembled here
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##INFO=<ID=SamlaFilter,Number=.,Type=String,Description=\"Samla filtering and case\">";
    header_vcf->addHeaderLine(line.str());

    line.str("");
    line << "##INFO=<ID=SamlaCase,Number=.,Type=String,Description=\"Complete annotation of Samla case decisions\">";
    header_vcf->addHeaderLine(line.str());

    Headers.add(Header_FILTER, "PASS",                     true,  "Site passed Samla tests, could reliably be called as either variant or reference");
    Headers.add(Header_FILTER, "FAIL",                     true,  "Site failed Samla tests, could not reliably be called as either variant or reference");
    Headers.add(Header_FILTER, "INFER_UNCALLED_VARIANT",   true,  "Samla method 'gwa', site passed quality threshold for variant but no variant available, appears to be uncalled variant");
    if (! opt_gwa_vqsr_vqsr_normal)
        Headers.add(Header_FILTER, "VQSR_VQSR_Fail",       true,  "Both G and W variants are VQSR-filtered, automatic failure");
    if (! opt_gwa_lowqual_lowqual_normal)
        Headers.add(Header_FILTER, "LowQual_LowQual_Fail", true,  "Both G and W variants are LowQual-filtered, automatic failure");
    if (! opt_gwa_mixed_normal)
        Headers.add(Header_FILTER, "Mixed_Fail",           true,  "One of G and W variants are VQSR-filtered, other is LowQual-filtered, automatic failure");
    if (opt_gwa_force_consistency)
        Headers.add(Header_FILTER, "Inconsistent",         true,  "G, W and A do not agree as to variant/no-variant in potentially ambiguous cases");
    Headers.add(Header_FILTER, "culpritFail",              true,  "Combination of VQSR culprits incompatible for variant call");
    Headers.add(Header_FILTER, "culpritPass",              false, "Combination of VQSR culprits compatible for variant call; not a filter, output on option");
    Headers.add(Header_FILTER, "quality:",                 false, "Samla method 'gwa' set of quality and contextual quality values for G, W and A; not a filter, output on option");
    Headers.add(Header_FILTER, "culprits:",                false, "Samla method 'gwa' set of VQSR culprits for G, W and A; not a filter, output on option");

    Headers.add(Header_FILTER, "noA",      false, "Samla method 'gwa', no variant in all-reads VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "snpA",     false, "Samla method 'gwa', variant in all-reads VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "noG",      false, "Samla method 'gwa', no variant in genomic VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "snpG",     false, "Samla method 'gwa', variant in genomic VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "noW",      false, "Samla method 'gwa', no variant in WGA VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "snpW",     false, "Samla method 'gwa', variant in WGA VCF; not a filter, output on option");
    Headers.add(Header_FILTER, "vqsrG",    false, "Samla method 'gwa', variant in genomic VCF VQSR-filtered; not a filter, output on option");
    Headers.add(Header_FILTER, "lowqualG", false, "Samla method 'gwa', variant in genomic VCF is LowQual; not a filter, output on option");
    Headers.add(Header_FILTER, "vqsrW",    false, "Samla method 'gwa', variant in WGA VCF VQSR-filtered; not a filter, output on option");
    Headers.add(Header_FILTER, "lowqualW", false, "Samla method 'gwa', variant in WGA VCF is LowQual; not a filter, output on option");

    ss.str(""); ss << "noQ" << opt_gwa_lowqual_quality;
    Headers.add(Header_FILTER, ss.str(),   true,  "Summed direct qualities of variants for site failed quality threshold");
    ss.str(""); ss << "Q" << opt_gwa_lowqual_quality;
    Headers.add(Header_FILTER, ss.str(),   false, "Summed direct qualities of variants for site passed quality threshold; not a filter, output on option");
    ss.str(""); ss << "noCQ" << opt_gwa_lowqual_quality;
    Headers.add(Header_FILTER, ss.str(),   true,  "Summed contextual qualities of variants for site failed quality threshold");
    ss.str(""); ss << "CQ" << opt_gwa_lowqual_quality;
    Headers.add(Header_FILTER, ss.str(),   false, "Summed contextual qualities of variants for site passed quality threshold; not a filter, output on option");
    ss.str(""); ss << "noQR" << opt_gwa_lowqual_quality_ref;
    Headers.add(Header_FILTER, ss.str(),   true,  "Summed qualities at site failed quality threshold to call as reference");
    ss.str(""); ss << "QR" << opt_gwa_lowqual_quality_ref;
    Headers.add(Header_FILTER, ss.str(),   false, "Summed qualities at site passed quality threshold to call as reference; not a filter, output on option");


    Headers.add(Header_FILTER, "GWA1_A",   false, "Samla method 'gwa' case 1, variant output is derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA1_G",   false, "Samla method 'gwa' case 1, variant output is derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA1_W",   false, "Samla method 'gwa' case 1, variant output is derived from W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA2_G_G", false, "Samla method 'gwa' case 2, variant called in G, output is derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA2_G_A", false, "Samla method 'gwa' case 2, variant called in G, output is derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA2_W_W", false, "Samla method 'gwa' case 2, variant called in W, output is derived from W; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA2_W_A", false, "Samla method 'gwa' case 2, variant called in W, output is derived from A; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA3_G",   false, "Samla method 'gwa' case 3, output is derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA3_W",   false, "Samla method 'gwa' case 3, output is derived from W; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA3_A",   false, "Samla method 'gwa' case 3, output is derived from A; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA4_A",   false, "Samla method 'gwa' case 4, output is derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA4_G",   false, "Samla method 'gwa' case 4, output is derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA4_W",   false, "Samla method 'gwa' case 4, output is derived from W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA5_A",   false, "Samla method 'gwa' case 5, output is derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA5_G",   false, "Samla method 'gwa' case 5, output is derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA5_W",   false, "Samla method 'gwa' case 5, output is derived from W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA6_G_A", false, "Samla method 'gwa' case 6, no variant in A and G, VQSR'd variant in W, output derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA6_G_G", false, "Samla method 'gwa' case 6, no variant in G, VQSR'd variant in W, variant in A, output derived from higher-quality G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA6_G_W", false, "Samla method 'gwa' case 6, no variant in G, VQSR'd variant in W, variant in A, output derived from higher-quality W; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA6_W_A", false, "Samla method 'gwa' case 6, no variant in A and W, VQSR'd variant in G, output derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA6_W_G", false, "Samla method 'gwa' case 6, no variant in W, VQSR'd variant in G, variant in A, output derived from higher-quality G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA6_W_W", false, "Samla method 'gwa' case 6, no variant in W, VQSR'd variant in G, variant in A, output derived from higher-quality W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA7_G_A", false, "Samla method 'gwa' case 7, no variant in A and G, LowQual variant in W, output derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA7_G_G", false, "Samla method 'gwa' case 7, no variant in G, LowQual variant in W, variant in A, output derived from higher-quality G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA7_G_W", false, "Samla method 'gwa' case 7, no variant in G, LowQual variant in W, variant in A, output derived from higher-quality W; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA7_W_A", false, "Samla method 'gwa' case 7, no variant in A and W, LowQual variant in G, output derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA7_W_G", false, "Samla method 'gwa' case 7, no variant in W, LowQual variant in G, variant in A, output derived from higher-quality G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA7_W_W", false, "Samla method 'gwa' case 7, no variant in W, LowQual variant in G, variant in A, output derived from higher-quality W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA8_G",   false, "Samla method 'gwa' case 8, no variant in W, variant in G, output derived from G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA8_W",   false, "Samla method 'gwa' case 8, no variant in G, variant in W, output derived from W; not a filter, output on option");

    Headers.add(Header_FILTER, "GWA9_A",   false, "Samla method 'gwa' case 9, no variant in A, G and W, output derived from A; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA9_G",   false, "Samla method 'gwa' case 9, no variant in G and W, variant in A, output derived from higher-quality G; not a filter, output on option");
    Headers.add(Header_FILTER, "GWA9_W",   false, "Samla method 'gwa' case 9, no variant in G and W, variant in A, output derived from higher-quality W; not a filter, output on option");

    Headers.fill_vcf_header(header_vcf, opt_full_filter_annotate);
}


int 
main(int argc, char* argv[]) {

    VcfStripmine vcfmine;

    processCommandLine(argc, argv, vcfmine);

    // determine which VCF to use for header lines and augment it
    VariantCallFile * header_vcf;
    if (samla_method == "gwa") {
        header_vcf = vcfmine.get_VariantCallFilePtr(vcf_all);
    } else if (samla_method == "default") {
        header_vcf = vcfmine.get_VariantCallFilePtr(vcfmine.first_vcf);
    } else {
        header_vcf = vcfmine.get_VariantCallFilePtr(vcfmine.first_vcf);
        cerr << "don't know the source of the initial VCF header, defaulting to the first" << endl;
    }
    prepareVcfHeader(header_vcf);

    cout << header_vcf->header << endl;  // output the header

    // initiate the stripmine, which fetches the first variant from each file
    if (! vcfmine.initiate()) 
        exitusage("Could not initiate stripmine");

    VcfStripmine::VariantConstPtr_vector vars;

    while (vcfmine.get(vars)) {

        // for each common set of positions containing variants in one or more VCF

        Variant v_ANS;

        if (samla_method == "gwa") {

            if (! method_gwa(v_ANS, vars))
                exitmessage("failure during method_gwa()");

            cout << v_ANS << endl;  // output the variant

        } else {

            exitmessage("unknown method '", samla_method, "'");

        }

        ++num_variants;
        if (opt_progress && (num_variants % opt_progress == 0)) {
            cerr << NAME "::main(): processed " << num_variants << " positions, last was " 
                << vars[0]->sequenceName << ":" << vars[0]->position << endl;
        }

    }

    return(0);
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
    ss << "quality:";
    ss << "G:" << v_Gen.quality << "/" << showpos << (v_Gen.quality - qualwindow_Gen.mean()) << noshowpos << ":";
    ss << "W:" << v_Wga.quality << "/" << showpos << (v_Wga.quality - qualwindow_Wga.mean()) << noshowpos << ":";
    ss << "A:" << v_All.quality << "/" << showpos << (v_All.quality - qualwindow_All.mean()) << noshowpos;
    return(ss.str());
}


enum culpritCase_t { culprit_fail = 0, culprit_none, culprit_pass_A, culprit_pass_G = culprit_pass_A, culprit_pass_W, culprit_pass_GW };

static culpritCase_t
culpritCase(Variant& v_Gen, Variant& v_Wga) {
    //if (! v_Gen.info.count("culprit") || !v_Wga.info.count("culprit")) {
    //    cerr << "culpritCase(Gen,Wga): no culprit when we expected one" << endl; exit(1);
    //}
    string g = v_Gen.info.count("culprit") ? v_Gen.info["culprit"][0] : "none";
    string w = v_Wga.info.count("culprit") ? v_Wga.info["culprit"][0] : "none";
    if ((g == "MQRankSum" || g == "ReadPosRankSum") && (w == "MQRankSum" || w == "ReadPosRankSum")) {
        return culprit_fail;
    } else if ((g == "DP" || g == "FS" || g == "QD") && (w == "MQRankSum" || w == "ReadPosRankSum")) {
        return culprit_pass_G;
    } else if ((g == "MQRankSum" || g == "ReadPosRankSum") && (w == "DP" || w == "FS" || w == "QD")) {
        return culprit_pass_W;
    } else if ((g == "DP" || g == "FS" || g == "QD") && (w == "DP" || w == "FS" || w == "QD")) {
        return culprit_pass_GW;
    } else {
        // note: "none" cases will fall here too
        cerr << "culpritCase(G,W): unhandled case g = " << g << " w = " << w << endl;
        exit(1);
    }
}

static culpritCase_t
culpritCase(Variant& v_A) {
    //if (v_A.info.count("culprit")) {
    //    cerr << "culpritCase(A): no culprit when we expected one" << endl; exit(1);
    //}
    string a = v_A.info.count("culprit") ? v_A.info["culprit"][0] : "none";
    if (a == "MQRankSum" || a == "ReadPosRankSum") {
        return culprit_fail;
    } else if (a == "DP" || a == "FS" || a == "QD") {
        return culprit_pass_A;
    } else {
        // note: "none" cases will fall here too
        cerr << "culpritCase(A): unhandled case a = " << a << endl;
        exit(1);
    }
}

string
generate_culprit_string(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    stringstream ss;
    ss << "culprit:";
    ss << "G:" << (v_Gen.info.count("culprit") ? v_Gen.info["culprit"][0] : "none") << ":";
    ss << "W:" << (v_Wga.info.count("culprit") ? v_Wga.info["culprit"][0] : "none") << ":";
    ss << "A:" << (v_All.info.count("culprit") ? v_All.info["culprit"][0] : "none");
    return(ss.str());
}

Variant
create_return_variant(Variant& v) {
    // clear any fields we need to
    Variant ANS = v;
    ANS.filter = "";
    ANS.info.erase("culprit");
    ANS.info.erase("VQSLOD");
    ANS.infoFlags.erase("NEGATIVE_TRAIN_SITE");
    ANS.infoFlags.erase("POSITIVE_TRAIN_SITE");
    return(ANS);
}

bool
is_GATK_variant(Variant& v) {
    // appropriate for "new" files from GATK, not for our modified alleles which will have PASS/FAIL for all sites
    return(v.filter == "PASS" || (v.filter == "." && v.alleles[1] != "."));
}

bool
is_GATK_variant_loose(Variant& v) {
    // appropriate for "new" files from GATK, not for our modified alleles which will have PASS/FAIL for all sites
    // differs from above since the variant might be filtered
    return(v.alleles[1] != ".");
}

bool
not_GATK_variant(Variant& v) {
    // appropriate for "new" files from GATK, not for our modified alleles which will have PASS/FAIL for all sites
    return((v.filter == "." || v.info["VariantType"][0] == "NO_VARIATION") && v.alleles[1] == ".");
}

bool
is_VQSR(Variant& v) {
    return(v.filter.substr(0, 4) == "VQSR");
}

bool
is_LowQual(Variant& v) {
    return(v.filter == "LowQual");
}

bool
is_VQSR_or_LowQual(Variant& v) {
    return(is_VQSR(v) || is_LowQual(v));
}

void
set_filter(Variant& v, const string& a) {
    if (! opt_filter_annotate && ! opt_full_filter_annotate && ! Headers.exists(a))
        exitmessage("*** set_filter() Unknown filter ", a, ", does not exist in Headers.");
    v.filter = a;
    v.info["SamlaFilter"].push_back(a);
}

void
set_first_filter(Variant& v, const string& a) {
    if (! opt_filter_annotate && ! opt_full_filter_annotate && ! Headers.exists(a))
        exitmessage("*** set_first_filter() Unknown filter ", a, ", does not exist in Headers.");
    if (v.filter == "" || v.filter == ".")
        v.filter = a;
    else
        v.filter = a + ";" + v.filter;
    v.info["SamlaFilter"].insert(v.info["SamlaFilter"].begin(), a);
}

void
annotate_filter(Variant& v, const string& a) {
    if (! opt_filter_annotate && ! opt_full_filter_annotate && ! Headers.exists(a))
        exitmessage("*** annotate_filter() Unknown filter ", a, ", does not exist in Headers.");
    if (opt_filter_annotate) {
        //v.addFilter(string("Samla") + a);
        v.addFilter(a);
    }
    v.info["SamlaFilter"].push_back(a);
}

void
annotate_case(Variant& v, const string& a) {
    v.info["SamlaCase"].push_back(a);
}

void
annotate_case_add_filter(Variant& v, const string& a) {
    annotate_filter(v, a);
    annotate_case(v, a);
}

void
annotate_case_add_full_filter(Variant& v, const string& a) {
    if (opt_full_filter_annotate)
        annotate_filter(v, a);
    annotate_case(v, a);
}

void
annotate_first_case(Variant& v, const string& a) {
    if (opt_full_filter_annotate)
        annotate_filter(v, a);
    v.info["SamlaCase"].insert(v.info["SamlaCase"].begin(), a);
}

static void
annotate_passed_q_filter(Variant& v, const string& tag, const double q) {
    stringstream ss;
    ss << tag << q;
    annotate_case_add_full_filter(v, ss.str());  // not included in FILTER field by default
}

static void
annotate_failed_q_filter(Variant& v, const string& tag, const double q) {
    stringstream ss;
    ss << tag << q;
    annotate_case_add_filter(v, ss.str());  // included in FILTER field
}

void
annotate_passed_Quality(Variant& v, const double q) { 
    annotate_passed_q_filter(v, "Q", q); 
}

void
annotate_failed_Quality(Variant& v, const double q) {
    annotate_failed_q_filter(v, "noQ", q);
}

void
annotate_passed_QualityRef(Variant& v, const double q) {
    annotate_passed_q_filter(v, "QR", q);
}

void
annotate_failed_QualityRef(Variant& v, const double q) {
    annotate_failed_q_filter(v, "noQR", q);
}

void
annotate_passed_ContextQuality(Variant& v, const double q) {
    annotate_passed_q_filter(v, "CQ", q);
}

void
annotate_failed_ContextQuality(Variant& v, const double q) {
    annotate_failed_q_filter(v, "noCQ", q);
}

// Implement gwa method, see comments below for details
bool method_gwa(Variant& v_ANS, VcfStripmine::VariantConstPtr_vector& vars) {

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
        if (skip_case && not_GATK_variant(v_Gen)
                      && not_GATK_variant(v_Wga)
                      && not_GATK_variant(v_All)) {  // no variants to consider
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

    // Variant v_ANS;

    if (is_GATK_variant(v_Gen) && is_GATK_variant(v_Wga)) {                 // CASE 1

        v_ANS = method_gwa_case1(v_Gen, v_Wga, v_All);

    } else if (is_GATK_variant(v_Gen) && is_VQSR_or_LowQual(v_Wga)) {       // CASE 2_G

        v_ANS = method_gwa_case2_G(v_Gen, v_Wga, v_All);

    } else if (is_VQSR_or_LowQual(v_Gen) && is_GATK_variant(v_Wga)) {       // CASE 2_W

        v_ANS = method_gwa_case2_W(v_Gen, v_Wga, v_All);

    } else if (is_VQSR(v_Gen) && is_VQSR(v_Wga)) {                          // CASE 3

        v_ANS = method_gwa_case3(v_Gen, v_Wga, v_All);

    } else if (is_LowQual(v_Gen) && is_LowQual(v_Wga)) {                    // CASE 4

        v_ANS = method_gwa_case4(v_Gen, v_Wga, v_All);

    } else if ((is_VQSR(v_Gen) && is_LowQual(v_Wga))
               ||
               (is_LowQual(v_Gen) && is_VQSR(v_Wga))) {                     // CASE 5

        v_ANS = method_gwa_case5(v_Gen, v_Wga, v_All);

    } else if (not_GATK_variant(v_Gen) && is_VQSR(v_Wga)) {                 // CASE 6_G

        v_ANS = method_gwa_case6_G(v_Gen, v_Wga, v_All);

    } else if (is_VQSR(v_Gen) && not_GATK_variant(v_Wga)) {                 // CASE 6_W

        v_ANS = method_gwa_case6_W(v_Gen, v_Wga, v_All);

    } else if (not_GATK_variant(v_Gen) && is_LowQual(v_Wga)) {              // CASE 7_G

        v_ANS = method_gwa_case7_G(v_Gen, v_Wga, v_All);

    } else if (is_LowQual(v_Gen) && not_GATK_variant(v_Wga)) {              // CASE 7_W

        v_ANS = method_gwa_case7_W(v_Gen, v_Wga, v_All);

    } else if (is_GATK_variant(v_Gen) && not_GATK_variant(v_Wga)) {         // CASE 8_G

        v_ANS = method_gwa_case8_G(v_Gen, v_Wga, v_All);

    } else if (not_GATK_variant(v_Gen) && is_GATK_variant(v_Wga)) {         // CASE 8_W

        v_ANS = method_gwa_case8_W(v_Gen, v_Wga, v_All);

    } else if (not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga) 
            && !skip_case) {                                                // CASE 9

        v_ANS = method_gwa_case9(v_Gen, v_Wga, v_All);

    } else {
        cerr << "=================  Unhandled case in method_gwa()" << endl;
        cerr << "**G " << v_Gen << endl << "**W " << v_Wga << endl << "**A " << v_All << endl;
        cerr << "---" << endl;
        exit(1);
    }

    string q = generate_gwa_qual_string(v_Gen, v_Wga, v_All);
    annotate_case_add_full_filter(v_ANS, q);

    // this is now done by create_return_variant()
    // clear INFO fields we shouldn't have
    // v_ANS.info.erase("culprit");

    // cout << v_ANS << endl;

    return(true);
}


// Note that VQSR can remove a SNP that otherwise has incredibly strong support, say for
// reasons of strand bias or mapping quality or read position bias.  Some of these reasons
// are expected to be stochastic and thus specific to a library, whereas others may be
// site-specific.
//
// FS: strand bias, library-specific
// DP: depth, library-specific
// QD: quality depth, library-specific
// MQRankSum: mapping quality rank sum, may be site-specific
// ReadPosRankSum: segment of read supporting variant, may be site-specific

/*********************** Current set of case labels

Case 1: Unambiguous variant for G and W.  What is emitted depends on whether A is
also a variant, and if not which of G and W has the highest quality.  Suffix indicates
which variant was emitted as the variant (_A, _G, _W).

GWA1_A : G, W, A are variants : Emit A as variant. PASS

    The unambiguous variant A is the emitted variant.  Quality is the summed qualities
    of G and W.

GWA1_G : G, W are variants, A is *not* a variant : Emit G as variant. PASS

    Since A was not a variant, we choose one of the others to emit.  Quality G
    >= W so G is the emitted variant.  Quality is quality of G.  If quality of G
    does not meet --gwa-quality, then site is a FAIL
    
GWA1_W : G, W are variants, A is *not* a variant : Emit W as variant. PASS

    Since A was not a variant, we choose one of the others to emit.  Quality G
    < W so W is the emitted variant.  Quality is quality of W.  If quality of W
    does not meet --gwa-quality, then site is a FAIL
    
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
        v_ANS = create_return_variant(v_All);
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        annotate_first_case(v_ANS, "GWA1_A");
        annotate_case(v_ANS, "Qual_G+W");
        annotate_case_add_full_filter(v_ANS, "snpA");
    } else {
        if (v_Gen.quality >= v_Wga.quality) { // favour v_Gen
            v_ANS = create_return_variant(v_Gen);  // using quality from v_Gen
            annotate_first_case(v_ANS, "GWA1_G");
            annotate_case(v_ANS, "Qual_G");
            annotate_case_add_full_filter(v_ANS, "snpG");
        } else {  // v_Gen.quality < v_Wga.quality, so favour v_Wga
            v_ANS = create_return_variant(v_Wga);  // using quality from v_Wga
            annotate_first_case(v_ANS, "GWA1_W");
            annotate_case(v_ANS, "Qual_W");
            annotate_case_add_full_filter(v_ANS, "snpW");
        }
        annotate_case_add_full_filter(v_ANS, "noA");
    }
    if (v_ANS.quality >= opt_gwa_quality) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_Quality(v_ANS, opt_gwa_quality);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_Quality(v_ANS, opt_gwa_quality);
    }
    return(v_ANS);
}


// -------- 2_G  Unambiguous SNP call for v_Gen, filtered for v_Wga


Variant
method_gwa_case2_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {

    // Combine strength from both if the filtered SNP is a variant, using
    // straight quality for the variant and contextual quality for the filtered
    // one if LowQual and a variant.  We don't know if we PASS until we have
    // all the qualities.
    if (DEBUG(2)) cout << "*** case 2_G method_gwa_case2_G()" << endl;
    bool W_LQ_not_variant = (is_LowQual(v_Wga) && not_GATK_variant(v_Wga));
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION" || W_LQ_not_variant) {
        v_ANS = create_return_variant(v_Gen);
        annotate_first_case(v_ANS, "GWA2_G_G");
        annotate_case_add_full_filter(v_ANS, "snpG");
        annotate_case_add_full_filter(v_ANS, "noA");
        if (DEBUG(2)) {
            cerr << "**2_G* v_Gen is PASS with other filtered but v_All is not a variant" << endl;
            cerr << "**2_G* G " << v_Gen << endl << "**2_G* W " << v_Wga << endl << "**2_G* A " << v_All << endl;
        }
    } else {
        v_ANS = create_return_variant(v_All);
        annotate_first_case(v_ANS, "GWA2_G_A");
        annotate_case_add_full_filter(v_ANS, "snpG");
        annotate_case_add_full_filter(v_ANS, "snpA");
    }
    if (is_VQSR(v_Wga) || ! W_LQ_not_variant) {
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        annotate_case(v_ANS, "Qual_G+W");
        annotate_case_add_full_filter(v_ANS, "snpW");
        if (is_VQSR(v_Wga)) {
            annotate_case_add_full_filter(v_ANS, "vqsrW");
        } else {
            annotate_case_add_full_filter(v_ANS, "lowqualW");
        }
    } else if (W_LQ_not_variant) { // LowQual and call to reference
        if (opt_gwa_enable_context_quality) {
            v_ANS.quality = v_Gen.quality + (v_Wga.quality - qualwindow_Wga.mean());
            annotate_case(v_ANS, "Qual_G+WContext");
        } else {
            annotate_case(v_ANS, "Qual_G");
        }
        annotate_case_add_full_filter(v_ANS, "noW");
        annotate_case_add_full_filter(v_ANS, "lowqualW");
    } else {
        cerr << "**2_G* unhandled v_Wga.filter case" << endl; exit(1);
    }
    if (v_ANS.quality >= opt_gwa_quality) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_Quality(v_ANS, opt_gwa_quality);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_Quality(v_ANS, opt_gwa_quality);
    }
    return(v_ANS);
}


// -------- 2_W  Unambiguous SNP call for v_Wga, filtered for v_Gen


Variant
method_gwa_case2_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {

    // Combine strength from both if the filtered SNP is a variant, using
    // straight quality for the variant and contextual quality for the filtered
    // one if LowQual and a variant.  We don't know if we PASS until we have
    // assembled all the qualities.
    if (DEBUG(2)) cout << "*** case 2_W method_gwa_case2_W()" << endl;
    bool G_LQ_not_variant = (is_LowQual(v_Gen) && not_GATK_variant(v_Gen));
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION" || G_LQ_not_variant) {
        v_ANS = create_return_variant(v_Wga);
        annotate_first_case(v_ANS, "GWA2_W_W");
        annotate_case_add_full_filter(v_ANS, "snpW");
        annotate_case_add_full_filter(v_ANS, "noA");
        if (DEBUG(2)) {
            cerr << "**2_W* v_Wga is PASS with other filtered but v_All is not a variant" << endl;
            cerr << "**2_W* G " << v_Gen << endl << "**2_G* W " << v_Wga << endl << "**2_G* A " << v_All << endl;
        }
    } else {
        v_ANS = create_return_variant(v_All);
        annotate_first_case(v_ANS, "GWA2_W_A");
        annotate_case_add_full_filter(v_ANS, "snpW");
        annotate_case_add_full_filter(v_ANS, "snpA");
    }
    if (is_VQSR(v_Gen) || ! G_LQ_not_variant) {
        v_ANS.quality = v_Wga.quality + v_Gen.quality;
        annotate_case(v_ANS, "Qual_W+G");
        annotate_case_add_full_filter(v_ANS, "snpG");
        if (is_VQSR(v_Gen)) {
            annotate_case_add_full_filter(v_ANS, "vqsrG");
        } else {
            annotate_case_add_full_filter(v_ANS, "lowqualG");
        }
    } else if (G_LQ_not_variant) { // LowQual and call to reference
        if (opt_gwa_enable_context_quality) {
            v_ANS.quality = v_Wga.quality + (v_Gen.quality - qualwindow_Gen.mean());
            annotate_case(v_ANS, "Qual_W+GContext");
        } else {
            annotate_case(v_ANS, "Qual_W");
        }
        annotate_case_add_full_filter(v_ANS, "noG");
        annotate_case_add_full_filter(v_ANS, "lowqualG");
    } else {
        cerr << "**2_W* unhandled v_Gen.filter case" << endl; exit(1);
    }
    if (v_ANS.quality >= opt_gwa_quality) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_Quality(v_ANS, opt_gwa_quality);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_Quality(v_ANS, opt_gwa_quality);
    }
    return(v_ANS);
}


// -------- 3  VQSR genotypes for both v_Gen and v_Wga, so we likely have a variant here


Variant
method_gwa_case3(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {

    if (DEBUG(2)) cout << "*** case 3 method_gwa_case3()" << endl;
    Variant v_ANS;

    // check to see if the VQSR reasons were "stochastic" (FS, QD, DP) and allow us to emit a variant
    culpritCase_t vqsr_culprits = culpritCase(v_Gen, v_Wga);

    if (opt_gwa_vqsr_vqsr_normal == false) {

        // fail all these cases
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA3_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "VQSR_VQSR_Fail");

    } else if (vqsr_culprits != culprit_fail) {

        if (v_All.info["VariantType"][0] == "NO_VARIATION") {

            // unfortunately A does not think we have a variant, pick the best of the others
            if (v_Gen.quality >= v_Wga.quality) {
                v_ANS = create_return_variant(v_Gen);
                annotate_first_case(v_ANS, "GWA3_G");
                annotate_case(v_ANS, "Qual_G");
            } else {
                v_ANS = create_return_variant(v_Wga);
                annotate_first_case(v_ANS, "GWA3_W");
                annotate_case(v_ANS, "Qual_W");
            }
            annotate_case_add_full_filter(v_ANS, "noA");

        } else if (v_All.info["VariantType"][0] == "SNP" || v_All.info["VariantType"][0] == "MULTIALLELIC_SNP") {

            // A contains a variant, update its quality to reflect the combination
            v_ANS = create_return_variant(v_All);
            v_ANS.quality = v_Gen.quality + v_Wga.quality;
            annotate_first_case(v_ANS, "GWA3_A");
            annotate_case(v_ANS, "Qual_G+W");
            annotate_case_add_full_filter(v_ANS, "snpA");

        } else {

            cerr << "**3* unhandled v_All VariantType case: " << v_All.info["VariantType"][0] << endl; 
            exit(1);

        }
        annotate_case_add_full_filter(v_ANS, "culpritPass");
        if (v_ANS.quality >= opt_gwa_vqsr_quality) {
            set_first_filter(v_ANS, "PASS");
            annotate_passed_Quality(v_ANS, opt_gwa_vqsr_quality);
        } else {
            set_first_filter(v_ANS, "FAIL");
            annotate_failed_Quality(v_ANS, opt_gwa_vqsr_quality);
        }

    } else { // culprit_fail

        // Accept whatever A contains, but fail the site.
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA3_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "culpritFail");

    }

    // annotate with the culprit values
    string q = generate_culprit_string(v_Gen, v_Wga, v_All);
    annotate_case_add_full_filter(v_ANS, q);

    return(v_ANS);
}


// -------- 4  LowQual genotypes for both v_Gen and v_Wga, so do we pass a threshold to emit a genotype?
//
// This should not be combined with case 3 because of its use of contextual quality, which case 3 does not use


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
    double q_Gen = (opt_gwa_enable_context_quality) ? abs(qdelta_Gen) : v_Gen.quality;
    double q_Wga = (opt_gwa_enable_context_quality) ? abs(qdelta_Wga) : v_Wga.quality;
    double q_site = q_Gen + q_Wga;
    bool   passedQuality = (q_site >= opt_gwa_lowqual_quality);
    bool   passedQualityRef = (v_Gen.quality + v_Wga.quality >= opt_gwa_lowqual_quality_ref);
    bool   all_not_variant = (not_GATK_variant(v_All) && not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga));
    bool   all_is_variant = (is_GATK_variant_loose(v_All) && is_GATK_variant_loose(v_Gen) && is_GATK_variant_loose(v_Wga));
    bool   G_W_not_variant = (not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga));
    bool   G_W_is_variant = (is_GATK_variant_loose(v_Gen) && is_GATK_variant_loose(v_Wga));
    if (DEBUG(1)) {
        fprintf(stderr, "4 G quality: %9.4f  lookback: %9.4f  delta: %9.4f  q_Gen: %9.4f\n", v_Gen.quality, qualwindow_Gen.mean(), qdelta_Gen, q_Gen);
        fprintf(stderr, "4 W quality: %9.4f  lookback: %9.4f  delta: %9.4f  q_Wga: %9.4f\n", v_Wga.quality, qualwindow_Wga.mean(), qdelta_Wga, q_Wga);
        fprintf(stderr, "4 A quality: %9.4f  lookback: %9.4f  delta: %9.4f\n", v_All.quality, qualwindow_All.mean(), qdelta_All);
        fprintf(stderr, "4                   q_site G + W delta: %9.4f ", q_site);
        cerr << "passedQuality " << (passedQuality ? "PASSED" : "FAILED") << " gwa-lowqual-quality=" << opt_gwa_lowqual_quality << endl;
        cerr << "passedQualityRef " << (passedQualityRef ? "PASSED" : "FAILED") << " gwa-lowqual-quality-ref=" << opt_gwa_lowqual_quality_ref << endl;
    }

    if (opt_gwa_lowqual_lowqual_normal == false) {

        //fail all these cases
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA4_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "LowQual_LowQual_Fail");

    } else if (opt_gwa_force_consistency && ! (all_not_variant || all_is_variant)) {

        // fail all these cases
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA4_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "Inconsistent");

    } else if (all_not_variant && passedQualityRef) {

#define annotate_4_passedQualityRef(__S__) \
            { \
                annotate_first_case(v_ANS, __S__); \
                if (opt_gwa_enable_context_quality) { \
                    annotate_failed_ContextQuality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_GContext+WContext"); \
                } else { \
                    annotate_failed_Quality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_G+W"); \
                } \
                annotate_passed_QualityRef(v_ANS, opt_gwa_lowqual_quality_ref); \
            }

        // we have no variant anywhere and passed QualityRef

        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "PASS");
        v_ANS.quality = q_site;
        annotate_4_passedQualityRef("GWA4_A");
        annotate_case_add_full_filter(v_ANS, "noG");
        annotate_case_add_full_filter(v_ANS, "noW");
        annotate_case_add_full_filter(v_ANS, "noA");

    } else if (passedQuality && (G_W_not_variant || G_W_is_variant)) {

#define annotate_4_passedQuality(__S__) \
            { \
                annotate_first_case(v_ANS, __S__); \
                if (opt_gwa_enable_context_quality) { \
                    annotate_passed_ContextQuality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_GContext+WContext"); \
                } else { \
                    annotate_passed_Quality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_G+W"); \
                } \
            }

        // There is a variant here, do we have one available in A?

        if (v_All.info["VariantType"][0] == "SNP" || v_All.info["VariantType"][0] == "MULTIALLELIC_SNP") {

            // use the variant in A
            v_ANS = create_return_variant(v_All);
            set_filter(v_ANS, "PASS");
            v_ANS.quality = q_site;
            annotate_4_passedQuality("GWA4_A");
            annotate_case_add_full_filter(v_ANS, "snpA");

        } else if (v_All.info["VariantType"][0] == "NO_VARIATION") {

            // if G or W contain a variant, use the one of higher quality
            // if neither contains a variant, we have to fail the site

            if (opt_gwa_enable_context_quality && not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga)) {

                // no variant is available but we think there should be, fail the site

                v_ANS = create_return_variant(v_All);
                set_filter(v_ANS, "FAIL");
                v_ANS.quality = q_site;
                annotate_4_passedQuality("GWA4_A");
                annotate_case_add_filter(v_ANS, "INFER_UNCALLED_VARIANT");
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "noW");
                annotate_case_add_full_filter(v_ANS, "noA");

            } else {

                if (G_W_is_variant) {

                    // both are available, pick one of higher quality and make its quality the contextual quality
                    if (v_Gen.quality >= v_Wga.quality) {
                        v_ANS = create_return_variant(v_Gen);
                        set_filter(v_ANS, "PASS");
                        annotate_4_passedQuality("GWA4_G");
                    } else {
                        v_ANS = create_return_variant(v_Wga);
                        set_filter(v_ANS, "PASS");
                        annotate_4_passedQuality("GWA4_W");
                    }
                    // v_ANS.quality = q_site;
                    annotate_case_add_full_filter(v_ANS, "snpG");
                    annotate_case_add_full_filter(v_ANS, "snpW");
                    annotate_case_add_full_filter(v_ANS, "noA");

                } else if (! not_GATK_variant(v_Gen)) {

                    // v_Gen has a variant
                    v_ANS = create_return_variant(v_Gen);
                    set_filter(v_ANS, "PASS");
                    // v_ANS.quality = q_site;
                    annotate_4_passedQuality("GWA4_G");
                    annotate_case_add_full_filter(v_ANS, "snpG");
                    annotate_case_add_full_filter(v_ANS, "noW");
                    annotate_case_add_full_filter(v_ANS, "noA");

                } else if (! not_GATK_variant(v_Wga)) {

                    // v_Wga has a variant
                    v_ANS = create_return_variant(v_Wga);
                    set_filter(v_ANS, "PASS");
                    // v_ANS.quality = q_site;
                    annotate_4_passedQuality("GWA4_W");
                    annotate_case_add_full_filter(v_ANS, "snpW");
                    annotate_case_add_full_filter(v_ANS, "noG");
                    annotate_case_add_full_filter(v_ANS, "noA");

                } else {

                    cerr << "**4* unhandled case within passedQuality" << endl; 
                    exit(1);

                }
            }

        } else {

            cerr << "**4* unhandled v_All VariantType case: " << v_All.info["VariantType"][0] << endl; 
            exit(1);

        }

    } else if (passedQualityRef) {

        // We do not think we should have a variant here, we think we should match the reference, can we do that?

        if (v_All.info["VariantType"][0] == "NO_VARIATION") {

            // use the no-variant in A, with quality is sum of context quality
            v_ANS = create_return_variant(v_All);
            set_filter(v_ANS, "PASS");
            v_ANS.quality = q_site;
            annotate_4_passedQualityRef("GWA4_A");
            annotate_case_add_full_filter(v_ANS, "noA");

        } else if (v_All.info["VariantType"][0] == "SNP" || v_All.info["VariantType"][0] == "MULTIALLELIC_SNP") {

            // if G or W do not contain a variant, use the one of higher quality
            // if both contains a variant, we have to fail the site

            if (G_W_not_variant) {

                // both non-variants are available, pick one of higher quality and make its quality the contextual quality
                if (v_Gen.quality >= v_Wga.quality) {
                    v_ANS = create_return_variant(v_Gen);
                    set_filter(v_ANS, "PASS");
                    v_ANS.quality = q_site;
                    annotate_4_passedQualityRef("GWA4_G");
                } else {
                    v_ANS = create_return_variant(v_Wga);
                    set_filter(v_ANS, "PASS");
                    v_ANS.quality = q_site;
                    annotate_4_passedQualityRef("GWA4_W");
                }
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "noW");
                annotate_case_add_full_filter(v_ANS, "snpA");

            } else if (not_GATK_variant(v_Gen)) {

                // v_Gen is non-variant, assuming that v_Wga is a (mistaken) variant
                v_ANS = create_return_variant(v_Gen);
                set_filter(v_ANS, "PASS");
                annotate_4_passedQualityRef("GWA4_G");
                annotate_case_add_full_filter(v_ANS, "snpW");
                annotate_case_add_full_filter(v_ANS, "noG");

            } else if (not_GATK_variant(v_Wga)) {

                // v_Wga is non-variant, assuming that v_Gen is a (mistaken) variant
                v_ANS = create_return_variant(v_Wga);
                set_filter(v_ANS, "PASS");
                annotate_4_passedQualityRef("GWA4_W");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "noW");

            } else if (G_W_is_variant) {

                // only variants available, it cannot be reference, fail the site
                v_ANS = create_return_variant(v_All);
                set_filter(v_ANS, "FAIL");
                v_ANS.quality = q_site;
                annotate_first_case(v_ANS, "GWA4_A");
                if (opt_gwa_enable_context_quality) {
                    annotate_failed_ContextQuality(v_ANS, opt_gwa_lowqual_quality);
                    annotate_case(v_ANS, "Qual_GContext+WContext");
                } else {
                    annotate_failed_Quality(v_ANS, opt_gwa_lowqual_quality);
                    annotate_case(v_ANS, "Qual_G+W");
                }
                annotate_case_add_full_filter(v_ANS, "snpA");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "snpW");


            } else {

                cerr << "**4* unhandled case within passedQualityRef" << endl; 
                exit(1);

            }

            // annotate_case(v_ANS, "IgnoreSnp_A");

        } else {

            cerr << "**4* unhandled v_All VariantType case: " << v_All.info["VariantType"][0] << endl; 
            exit(1);

        }

    } else {

        // We can't decide what to do, fail the site
        // Annotate re: variants but don't try to pick the right one, we don't know which that could be
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA4_A");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_All) ? "snpA" : "noA");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_Gen) ? "snpG" : "noG");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_Wga) ? "snpW" : "noW");
        v_ANS.quality = q_site;
        if (opt_gwa_enable_context_quality) {
            annotate_failed_ContextQuality(v_ANS, opt_gwa_lowqual_quality);
            annotate_case(v_ANS, "Qual_GContext+WContext");
        } else {
            annotate_failed_Quality(v_ANS, opt_gwa_lowqual_quality);
            annotate_case(v_ANS, "Qual_G+W");
        }
        annotate_failed_QualityRef(v_ANS, opt_gwa_lowqual_quality_ref);

    }

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
    double q_Gen = (v_Gen.filter == "LowQual" && opt_gwa_enable_context_quality) ? abs(qdelta_Gen) : v_Gen.quality;
    double q_Wga = (v_Wga.filter == "LowQual" && opt_gwa_enable_context_quality) ? abs(qdelta_Wga) : v_Wga.quality;
    culpritCase_t vqsr_culprit = (v_Gen.filter == "LowQual") ? culpritCase(v_Wga) : culpritCase(v_Gen);
    string culprit_string = generate_culprit_string(v_Gen, v_Wga, v_All);
    double q_site = q_Gen + q_Wga;
    bool   passedQuality = (q_site >= opt_gwa_mixed_quality);
    bool   passedQualityRef = (v_Gen.quality + v_Wga.quality >= opt_gwa_mixed_quality_ref);
    bool   all_not_variant = (not_GATK_variant(v_All) && not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga));
    bool   all_is_variant = (is_GATK_variant_loose(v_All) && is_GATK_variant_loose(v_Gen) && is_GATK_variant_loose(v_Wga));
    bool   G_W_not_variant = (not_GATK_variant(v_Gen) && not_GATK_variant(v_Wga));
    bool   G_W_is_variant = (is_GATK_variant_loose(v_Gen) && is_GATK_variant_loose(v_Wga));
    if (DEBUG(1)) {
        fprintf(stderr, "5 G quality: %9.4f  lookback: %9.4f  delta: %9.4f  q_Gen: %9.4f\n", v_Gen.quality, qualwindow_Gen.mean(), qdelta_Gen, q_Gen);
        fprintf(stderr, "5 W quality: %9.4f  lookback: %9.4f  delta: %9.4f  q_Wga: %9.4f\n", v_Wga.quality, qualwindow_Wga.mean(), qdelta_Wga, q_Wga);
        fprintf(stderr, "5 A quality: %9.4f  lookback: %9.4f  delta: %9.4f\n", v_All.quality, qualwindow_All.mean(), qdelta_All);
        fprintf(stderr, "5                   q_site G + W delta: %9.4f ", q_site);
        cerr << "passedQuality " << (passedQuality ? "PASSED" : "FAILED") << " gwa-mixed-quality=" << opt_gwa_mixed_quality << endl;
        cerr << "passedQualityRef " << (passedQualityRef ? "PASSED" : "FAILED") << " gwa-mixed-quality-ref=" << opt_gwa_mixed_quality_ref << endl;
        cerr << "vqsr_culprit = " << vqsr_culprit << "  culprit_string = " << culprit_string << endl;
    }

    if (opt_gwa_mixed_normal == false) {

        //fail all these cases
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA5_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "Mixed_Fail");

    } else if (opt_gwa_force_consistency && ! (all_not_variant || all_is_variant)) {

        // fail all these cases
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA5_A");
        annotate_case(v_ANS, "Qual_A");
        annotate_case_add_filter(v_ANS, "Inconsistent");

    } else if (vqsr_culprit != culprit_fail && passedQuality) {

        // There is a variant here, do we have one available in A?

#define annotate_5_passedQuality(__S__) \
            { \
                annotate_first_case(v_ANS, __S__); \
                if (opt_gwa_enable_context_quality) { \
                    annotate_passed_ContextQuality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_GContext+WContext"); \
                } else { \
                    annotate_passed_Quality(v_ANS, opt_gwa_lowqual_quality); \
                    annotate_case(v_ANS, "Qual_G+W"); \
                } \
                annotate_case_add_full_filter(v_ANS, "culpritPass"); \
                annotate_case_add_full_filter(v_ANS, culprit_string); \
            }

        // There is a variant here, do we have one available in A?

        if (v_All.info["VariantType"][0] == "SNP" || v_All.info["VariantType"][0] == "MULTIALLELIC_SNP") {

            // use the variant in A
            v_ANS = create_return_variant(v_All);
            set_filter(v_ANS, "PASS");
            v_ANS.quality = q_site;
            annotate_5_passedQuality("GWA5_A");
            annotate_case_add_full_filter(v_ANS, "snpA");

        } else if (v_All.info["VariantType"][0] == "NO_VARIATION") {

            // if G or W contain a variant, use the one of higher quality
            // if neither contains a variant, we have to fail the site

            if (G_W_is_variant) {

                // both are available, pick one of higher quality and adjust its quality
                if (v_Gen.quality >= v_Wga.quality) {
                    v_ANS = create_return_variant(v_Gen);
                    annotate_5_passedQuality("GWA5_G");
                } else {
                    v_ANS = create_return_variant(v_Wga);
                    annotate_5_passedQuality("GWA5_W");
                }
                set_first_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_case_add_full_filter(v_ANS, "noA");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "snpW");

            } else if (is_GATK_variant_loose(v_Gen)) {

                // v_Gen has a variant
                v_ANS = create_return_variant(v_Gen);
                set_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_5_passedQuality("GWA5_G");
                annotate_case_add_full_filter(v_ANS, "noA");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "noW");

            } else if (is_GATK_variant_loose(v_Wga)) {

                // v_Wga has a variant
                v_ANS = create_return_variant(v_Wga);
                set_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_5_passedQuality("GWA5_W");
                annotate_case_add_full_filter(v_ANS, "noA");
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "snpW");

            } else if (G_W_not_variant) {

                // no variant is available, fail the site
                v_ANS = create_return_variant(v_All);
                set_filter(v_ANS, "FAIL");
                v_ANS.quality = q_site;
                annotate_5_passedQuality("GWA5_A");
                annotate_case_add_filter(v_ANS, "INFER_UNCALLED_VARIANT");
                annotate_case_add_full_filter(v_ANS, "noA");
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "noW");

            } else {
                cerr << "**5* unhandled case within vqsr_culprit and passedQuality" << endl; 
                exit(1);
            }

        } else {
            cerr << "**5* unhandled v_All VariantType case: " << v_All.info["VariantType"][0] << endl; 
            exit(1);
        }

    } else if (vqsr_culprit != culprit_fail && passedQualityRef) { // we failed passedQuality or the culprit test, but we can match the ref

#define annotate_5_passedQualityRef(__S__) \
            { \
                annotate_first_case(v_ANS, __S__); \
                annotate_passed_QualityRef(v_ANS, opt_gwa_mixed_quality_ref); \
                if (opt_gwa_enable_context_quality) { \
                    annotate_failed_ContextQuality(v_ANS, opt_gwa_mixed_quality); \
                    annotate_case(v_ANS, "Qual_GContext+WContext"); \
                } else { \
                    annotate_failed_Quality(v_ANS, opt_gwa_mixed_quality); \
                    annotate_case(v_ANS, "Qual_G+W"); \
                } \
                if (vqsr_culprit == culprit_fail) { \
                    annotate_case_add_filter(v_ANS, "culpritFail"); \
                } \
                annotate_case_add_full_filter(v_ANS, culprit_string); \
            }

        // We do not think we should have a variant here, we think we should match the reference, can we do that?

        if (v_All.info["VariantType"][0] == "NO_VARIATION") {

            // use the no-variant in A, with quality is sum of context quality
            v_ANS = create_return_variant(v_All);
            set_filter(v_ANS, "PASS");
            v_ANS.quality = q_site;
            annotate_5_passedQualityRef("GWA5_A");
            annotate_case_add_full_filter(v_ANS, "noA");

        } else if (v_All.info["VariantType"][0] == "SNP" || v_All.info["VariantType"][0] == "MULTIALLELIC_SNP") {

            // if G or W do not contain a variant, use the one of higher quality
            // if both contains a variant, we have to fail the site

            if (G_W_not_variant) {
                // both non-variants are available, pick one of higher quality and make its quality the contextual quality
                if (v_Gen.quality >= v_Wga.quality) {
                    v_ANS = create_return_variant(v_Gen);
                    annotate_5_passedQualityRef("GWA5_G");
                } else {
                    v_ANS = create_return_variant(v_Wga);
                    annotate_5_passedQualityRef("GWA5_W");
                }
                set_first_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_case_add_full_filter(v_ANS, "snpA");
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "noW");

            } else if (not_GATK_variant(v_Gen)) {

                // v_Gen is non-variant, assuming that v_Wga is a (mistaken) variant
                v_ANS = create_return_variant(v_Gen);
                set_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_5_passedQualityRef("GWA5_G");
                annotate_case_add_full_filter(v_ANS, "snpA");
                annotate_case_add_full_filter(v_ANS, "noG");
                annotate_case_add_full_filter(v_ANS, "snpW");

            } else if (not_GATK_variant(v_Wga)) {

                // v_Wga is non-variant, assuming that v_Gen is a (mistaken) variant
                v_ANS = create_return_variant(v_Wga);
                set_filter(v_ANS, "PASS");
                v_ANS.quality = q_site;
                annotate_5_passedQualityRef("GWA5_W");
                annotate_case_add_full_filter(v_ANS, "snpA");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "noW");

            } else if (G_W_is_variant) {

                // only variants available, fail the site
                v_ANS = create_return_variant(v_All);
                set_filter(v_ANS, "FAIL");
                v_ANS.quality = q_site;
                annotate_5_passedQualityRef("GWA5_A");
                annotate_case_add_full_filter(v_ANS, "snpA");
                annotate_case_add_full_filter(v_ANS, "snpG");
                annotate_case_add_full_filter(v_ANS, "snpW");

            } else {
                cerr << "**5* unhandled case within passedQualityRef" << endl; 
                exit(1);
            }

        } else {
            cerr << "**5* unhandled v_All VariantType case: " << v_All.info["VariantType"][0] << endl; 
            exit(1);
        }

    } else {

        // We can't decide what to do, fail the site
        // Annotate re: variants but don't try to pick the right one, we don't know which that could be
        v_ANS = create_return_variant(v_All);
        set_filter(v_ANS, "FAIL");
        annotate_first_case(v_ANS, "GWA5_A");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_All) ? "snpA" : "noA");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_Gen) ? "snpG" : "noG");
        annotate_case_add_full_filter(v_ANS, is_GATK_variant_loose(v_Wga) ? "snpW" : "noW");
        v_ANS.quality = q_site;
        if (opt_gwa_enable_context_quality) {
            annotate_failed_ContextQuality(v_ANS, opt_gwa_mixed_quality);
            annotate_case(v_ANS, "Qual_GContext+WContext");
        } else {
            annotate_failed_Quality(v_ANS, opt_gwa_mixed_quality);
            annotate_case(v_ANS, "Qual_G+W");
        }
        annotate_failed_QualityRef(v_ANS, opt_gwa_mixed_quality_ref);
        if (vqsr_culprit == culprit_fail) {
            annotate_case_add_filter(v_ANS, "culpritFail");
        } else {
            annotate_case_add_full_filter(v_ANS, "culpritPass");
        }
        annotate_case_add_full_filter(v_ANS, culprit_string);

    }

    return(v_ANS);
}


// -------- 6_G  No variant in v_Gen, filtered variant (VQSR) in v_Wga.  Emit no variant, call ref from v_Gen


Variant
method_gwa_case6_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 6_G method_gwa_case6_G()" << endl;
    Variant v_ANS;
    v_ANS = create_return_variant(v_Gen);
    annotate_first_case(v_ANS, "GWA6_G");
    annotate_case(v_ANS, "Qual_G");
    annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_All) ? "noA" : "snpA"));
    annotate_case_add_full_filter(v_ANS, "noG");
    annotate_case_add_full_filter(v_ANS, "snpW");
    annotate_case_add_full_filter(v_ANS, "vqsrW");
    if (DEBUG(2)) {
        cerr << "**6_G** v_Gen is '.' with v_Wga LowQual but v_All *is* a variant" << endl;
        cerr << "**6_G** G " << v_Gen << endl << "**6_G** W " << v_Wga << endl << "**6_G** A " << v_All << endl;
    }
    if (v_ANS.quality >= opt_gwa_quality_ref) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_QualityRef(v_ANS, opt_gwa_quality_ref);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_QualityRef(v_ANS, opt_gwa_quality_ref);
    }
    return(v_ANS);
}


// -------- 6_W  Filtered variant (VQSR) in v_Gen, no variant in v_Wga.  Emit no variant.


Variant
method_gwa_case6_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 6_W method_gwa_case6_W()" << endl;
    Variant v_ANS;
    v_ANS = create_return_variant(v_Wga);
    annotate_first_case(v_ANS, "GWA6_W");
    annotate_case(v_ANS, "Qual_W");
    annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_All) ? "noA" : "snpA"));
    annotate_case_add_full_filter(v_ANS, "snpG");
    annotate_case_add_full_filter(v_ANS, "noW");
    annotate_case_add_full_filter(v_ANS, "vqsrG");
    if (DEBUG(2)) {
        cerr << "**6_W** v_Gen is LowQual with v_Wga '.' but v_All *is* a variant" << endl;
        cerr << "**6_W** G " << v_Gen << endl << "**6_W** W " << v_Wga << endl << "**6_W** A " << v_All << endl;
    }
    if (v_ANS.quality >= opt_gwa_quality_ref) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_QualityRef(v_ANS, opt_gwa_quality_ref);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_QualityRef(v_ANS, opt_gwa_quality_ref);
    }
    return(v_ANS);
}


// -------- 7_G  No variant in v_Gen, filtered site (LowQual) in v_Wga.  Emit no variant.


Variant
method_gwa_case7_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 7_G method_gwa_case7_G()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION" && not_GATK_variant(v_Wga)) {
        v_ANS = create_return_variant(v_All);
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        annotate_first_case(v_ANS, "GWA7_G_A");
        annotate_case(v_ANS, "Qual_G+W");
        annotate_case_add_full_filter(v_ANS, "noA");
        annotate_case_add_full_filter(v_ANS, "noG");
        annotate_case_add_full_filter(v_ANS, "noW");
    } else {
        // regardless of whatever reason we failed above, we only go with v_Gen and its quality here
        v_ANS = create_return_variant(v_Gen);
        annotate_first_case(v_ANS, "GWA7_G_G");
        annotate_case(v_ANS, "Qual_G");
        annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_All) ? "noA" : "snpA"));
        annotate_case_add_full_filter(v_ANS, "noG");
        annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_Wga) ? "noW" : "snpW"));
        if (DEBUG(2)) {
            cerr << "**7_G** v_Gen is '.' with v_Wga LowQual but v_All *is* a variant" << endl;
            cerr << "**7_G** G " << v_Gen << endl << "**7_G** W " << v_Wga << endl << "**7_G** A " << v_All << endl;
        }
    }
    if (v_ANS.quality >= opt_gwa_quality_ref) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_QualityRef(v_ANS, opt_gwa_quality_ref);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_QualityRef(v_ANS, opt_gwa_quality_ref);
    }
    return(v_ANS);
}


// -------- 7_W  Filtered site (LowQual) in v_Gen, no variant in v_Wga.  Emit no variant.


Variant
method_gwa_case7_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 7_W method_gwa_case7_W()" << endl;
    Variant v_ANS;
    if (v_All.info["VariantType"][0] == "NO_VARIATION" && not_GATK_variant(v_Gen)) {
        v_ANS = create_return_variant(v_All);
        v_ANS.quality = v_Gen.quality + v_Wga.quality,
        annotate_first_case(v_ANS, "GWA7_W_A");
        annotate_case(v_ANS, "Qual_G+W");
        annotate_case_add_full_filter(v_ANS, "noA");
        annotate_case_add_full_filter(v_ANS, "noG");
        annotate_case_add_full_filter(v_ANS, "noW");
    } else {
        // regardless of whatever reason we failed above, we only go with v_Wga and its quality here
        v_ANS = create_return_variant(v_Wga);
        annotate_first_case(v_ANS, "GWA7_W_W");
        annotate_case(v_ANS, "Qual_W");
        annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_All) ? "noA" : "snpA"));
        annotate_case_add_full_filter(v_ANS, (not_GATK_variant(v_Wga) ? "noG" : "snpG"));
        annotate_case_add_full_filter(v_ANS, "noW");
        if (DEBUG(2)) {
            cerr << "**7_W** v_Gen is LowQual with v_Wga '.' but v_All *is* a variant" << endl;
            cerr << "**7_W** G " << v_Gen << endl << "**7_W** W " << v_Wga << endl << "**7_W** A " << v_All << endl;
        }
    }
    if (v_ANS.quality >= opt_gwa_quality_ref) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_QualityRef(v_ANS, opt_gwa_quality_ref);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_QualityRef(v_ANS, opt_gwa_quality_ref);
    }
    return(v_ANS);
}


// -------- 8_G  Variant in v_Gen, no variant in v_Wga.  Emit v_Gen variant.


Variant
method_gwa_case8_G(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 8_G method_gwa_case8_G()" << endl;
    Variant v_ANS = create_return_variant(v_Gen);
    annotate_first_case(v_ANS, "GWA8_G");
    annotate_case_add_full_filter(v_ANS, "snpG");
    if (opt_gwa_enable_context_quality) {
        v_ANS.quality += (v_Wga.quality - qualwindow_Wga.mean());
        annotate_case(v_ANS, "Qual_G+WContext");
    } else {
        annotate_case(v_ANS, "Qual_G");
    }
    if (v_ANS.quality >= opt_gwa_quality) {
        set_first_filter(v_ANS, "PASS");
        if (opt_gwa_enable_context_quality) {
            annotate_passed_ContextQuality(v_ANS, opt_gwa_quality);
        } else {
            annotate_passed_Quality(v_ANS, opt_gwa_quality);
        }
    } else {
        set_first_filter(v_ANS, "FAIL");
        if (opt_gwa_enable_context_quality) {
            annotate_failed_ContextQuality(v_ANS, opt_gwa_quality);
        } else {
            annotate_failed_Quality(v_ANS, opt_gwa_quality);
        }
    }
    return(v_ANS);
}


// -------- 8_W  No-variant in v_Gen, variant in v_Wga.  Emit v_Wga variant.


Variant
method_gwa_case8_W(Variant& v_Gen, Variant& v_Wga, Variant& v_All) {
    if (DEBUG(2)) cout << "*** case 8_W method_gwa_case8_W()" << endl;
    Variant v_ANS = create_return_variant(v_Wga);
    annotate_first_case(v_ANS, "GWA8_W");
    annotate_case_add_full_filter(v_ANS, "snpW");
    if (opt_gwa_enable_context_quality) {
        v_ANS.quality += (v_Gen.quality - qualwindow_Gen.mean());
        annotate_case(v_ANS, "Qual_W+GContext");
    } else {
        annotate_case(v_ANS, "Qual_W");
    }
    if (v_ANS.quality >= opt_gwa_quality) {
        set_first_filter(v_ANS, "PASS");
        if (opt_gwa_enable_context_quality) {
            annotate_passed_ContextQuality(v_ANS, opt_gwa_quality);
        } else {
            annotate_passed_Quality(v_ANS, opt_gwa_quality);
        }
    } else {
        set_first_filter(v_ANS, "FAIL");
        if (opt_gwa_enable_context_quality) {
            annotate_failed_ContextQuality(v_ANS, opt_gwa_quality);
        } else {
            annotate_failed_Quality(v_ANS, opt_gwa_quality);
        }
    }
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
        v_ANS = create_return_variant(v_All);
        annotate_first_case(v_ANS, "GWA9_A");
        v_ANS.quality = v_Gen.quality + v_Wga.quality;
        annotate_case(v_ANS, "Qual_G+W");
        annotate_case_add_full_filter(v_ANS, "noA");
    } else {
        // Favour the higher-quality no-variant if v_All disagrees
        // We will probably very rarely end up here
        if (v_Gen.quality >= v_Wga.quality) {
            v_ANS = create_return_variant(v_Gen);
            annotate_first_case(v_ANS, "GWA9_G");
            annotate_case(v_ANS, "Qual_G");
        } else {
            v_ANS = create_return_variant(v_Wga);
            annotate_first_case(v_ANS, "GWA9_W");
            annotate_case(v_ANS, "Qual_W");
        }
        annotate_case_add_full_filter(v_ANS, "snpA");
        if (DEBUG(2)) {
            cerr << "**9** v_Gen and v_Wga both no-variant, but v_All *is* a variant" << endl;
            cerr << "**9** G " << v_Gen << endl << "**9** W " << v_Wga << endl << "**9** A " << v_All << endl;
        }
    }
    annotate_case_add_full_filter(v_ANS, "noG");
    annotate_case_add_full_filter(v_ANS, "noW");
    if (v_ANS.quality >= opt_gwa_quality_ref) {
        set_first_filter(v_ANS, "PASS");
        annotate_passed_QualityRef(v_ANS, opt_gwa_quality_ref);
    } else {
        set_first_filter(v_ANS, "FAIL");
        annotate_failed_QualityRef(v_ANS, opt_gwa_quality_ref);
    }
    return(v_ANS);
}

