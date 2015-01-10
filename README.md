Samla
=====

Collect variants from parallel VCF files, apply heuristics, and produce a consensus set of variants.

Samla's `--method gwa` and `--method gwa-ksp` have been developed to produce consensus variants from those produced by different types of sequencing runs of the same individual.  Specifically, input is three GATK all-sites VCF files based on three types of sequencing data: (1) reads from whole-genome-shotgun sequencing; (2) reads from sequencing following whole-genome amplification; and (3) both sets of reads combined.

Coverage of data types (1) may be insufficient for high-confidence genotyping, hence the creation of (2).  However, coverage can be highly variable in (2).  Because the error models for variant calling are different between (1) and (2), it is not correct to use genotypes from (3) without further investigation.

Samla produces a consensus all-sites VCF file after examining the combined strength of support for both reference-matching and variant sites using all three VCF files.  The consensus VCF may contain modest gains in variants, on the order of a few percentage points, but both reference and variant sites will be of better quality: those with consistent support will be strengthened, while those with inconsistent support will be weakened.

More information is available via `samla --help`.

~~~~
Usage:   samla [options] -r refnames.txt <in1.vcf> [ in2.vcf ... ]

Collect VCF results and produce a consensus list of variants.

NOTE: samla 0.1.1-1-g9f77aea-52 is under active development.

     --references FILE            file containing reference names in order [REQUIRED]
                                  should be in the same order as the VCF files
     --output FILE                output file name [default is stdout]
     --debug INT                  debug info level INT [0]
     --progress INT               print variants processed mod INT [0]

     --filter-annotate            annotate FILTER field with additional method-specific information
     --full-filter-annotate       annotate FILTER field with even more method-specific information
     --no-filter-annotate         do not annotate FILTER field, use only PASS and filters indicating failure described in the VCF header [set]

     --method METHOD              use combining method 'METHOD', only 'gwa' and 'gwa-ksp' are implemented

'gwa' method options:

     --gwa-window INT                 Lookback window size for mean quality, max 50[5]
     --gwa-quality FLOAT              Minimum quality when combining both VQSR and LowQual variants [30]
                                      Specifying this option will set all the quality values below to the given value.
     --gwa-quality-ref FLOAT          Minimum quality when combining variants and call matches reference [20]
                                      Specifying this option will also set --gwa-lowqual-quality-ref and --gwa-mixed-quality-ref.
     --gwa-vqsr-quality FLOAT         Minimum quality when combining VQSR variants [30]
     --gwa-lowqual-quality FLOAT      Minimum quality when combining LowQual variants and call does not match reference [30]
     --gwa-lowqual-quality-ref FLOAT  Minimum quality to meet when combining LowQual variants and call matches reference [20]
     --gwa-mixed-quality FLOAT        Minimum quality when combining VQSR with LowQual variants [30]
     --gwa-mixed-quality-ref FLOAT    Minimum quality to meet when combining VQSR with LowQual variants and call matches reference [20]

     --gwa-force-consistency          Require G, W and A to agree on variant/no-variant for potentially ambiguous cases 4 and 5 [set]
     --gwa-no-force-consistency       Do not require G, W and A to agree on variant/no-variant for potentially ambiguous cases 4 and 5

     --gwa-disable-context-quality    Disables usage of context quality, qualities instead compared directly [set]
     --gwa-enable-context-quality     Enables usage of context quality
                                      NOTE: context quality is very experimental and may be incorrect, use at you rown risk

     --gwa-vqsr-vqsr-fail             Mark as FAIL all cases having both genomic and WGA VQSR-filtered variants [set]
     --gwa-vqsr-vqsr-normal           PASS/FAIL cases having both genomic and WGA VQSR-filtered variants depending on culprits
     --gwa-lowqual-lowqual-fail       Mark as FAIL all cases having both genomic and WGA LowQual-filtered variants
     --gwa-lowqual-lowqual-normal     PASS/FAIL cases having both genomic and WGA LowQual-filtered variants depending on culprits [set]
     --gwa-mixed-fail                 Mark as FAIL all cases having both a VQSR-filtered and a LowQual-filtered variant [set]
     --gwa-mixed-normal               PASS/FAIL cases having both a VQSR-filtered and a LowQual-filtered variant
                                      --gwa-vqsr-lowqual-fail and --gwa-vqsr-lowqual-normal are synonyms for the --gwa-mixed-* options
     --gwa-case8-emit-all             For case 8 (one of G/W has a variant while the other does not), emit the A call for the site [set]
                                      When the quality of one library (often W) is poor, this option may help reduce bias at the
                                      probable cost of a slight increase in failing case 8 sites for samples with better libraries.
     --gwa-case8-no-emit-all          For case 8 (one of G/W has a variant while the other does not), emit the variant


The 'gwa-ksp' method sets the following options. Options appearing after this may make further changes to option values.

     --method gwa
     --gwa-quality 30
     --gwa-quality-ref 20
     --gwa-force-consistency
     --gwa-disable-context-quality
     --gwa-vqsr-vqsr-fail
     --gwa-mixed-fail
     --gwa-case8-emit-all


For methods 'gwa' and 'gwa-ksp', all VCF files must be specified using these options:

     --vcf-genomic FILE           VCF file containing genomic calls
     --vcf-wga FILE               VCF file containing whole-genome-amplified calls
     --vcf-all FILE               VCF file containing all (pooled) calls


     -h | -? | --help             help

Version:     samla 0.1.1-1-g9f77aea-52

Compiler:    g++ (MacPorts gcc48 4.8.3_4) 4.8.3
Build flags: -Ivcflib -Wall -D_FILE_OFFSET_BITS=64 -O0 -D_WITH_DEBUG -ggdb -g3 -fvar-tracking-assignments -fno-inline -fno-inline-small-functions -fno-eliminate-unused-debug-types
~~~~

