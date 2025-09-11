# scope:
#     shared alignment wrapper used for both haploid reference genome and genotype graph alignment
# actions:
#     apply fastp quality filtering, adapter trimming, and read merging (when using short, paired reads)
#     align both high accuracy short and all long read sequences
# expects:
#     [source $MODULES_DIR/genotype/set_genotype_vars.sh] for genotype
#     source $MODULES_DIR/genome/set_genome_vars.sh
#     source $MODULES_DIR/align/set_alignment_vars.sh
#     source $MODULES_DIR/source/set_read_file_vars.sh
#     source $MODULES_DIR/agFree/set_agFree_vars.sh
# optional for all read types:
#     $FORCE_ALIGNMENT  [default: don't overwrite ALIGNER_OUTPUT_FILE]
# input:
#     $READ_1_FILES
#     $READ_2_FILES for paired read platforms
#     if ONT, reads subjected to trimming by modules/ONT/trim.pl (NOT by Dorado)
# output:
#     NAME_BAM_FILE
# where bam:
#     SEQ and QUAL are set to * unless needed for parsing variants (SEQ for SVs, QUAL for SNVs)
#     includes cs:Z: tag
#     includes custom tags:
#       xf:i: flag to record which alignment, reads, and events have candidate mosaic variants
# and where QNAME is modified with appended trailing pre-alignment tags:
#     QNAME:channel:trim5:trim3:mergeLevel:nRead1:nRead2
#     where:
#         QNAME is the event identifier provided by the sequencing platform
#         channel:trim5:trim3 is:
#             added with final values by basecall/ONT/trim.pl for ONT reads
#             added with null values 0:0:0 by align/do/prepare_fastq.pl for all other platforms
#         channel is the ONT channel, i.e., individual nanopore; 0 for all other platforms
#         trim5 is the number of adapter bases trimmed from the sequence 5' end (beginning of read1)
#         trim3 is the number of adapter bases trimmed from the sequence 3' end
#             trim values are NOT the internal, i.e., 3' adapters trimmed during paired read merging
#         mergeLevel:read1:read2 is added by align/do/adjust_merge_tags.pl
#         mergeLevel is 2 for fastp-merged reads, 0 for single or unmerged reads
#             mergeLevel may later be adjusted to 1 if subjected to downstream alignment-guided merging
#         nRead1 is the number of read1 bases present in the final merged read
#         nRead2 is the number of read2 bases present in the final merged read

#------------------------------------------------------------------
# set the product alignment file; abort silently if exists and not forced
#------------------------------------------------------------------
if [[ "$FORCE_ALIGNMENT" != "" && "$FORCE_ALIGNMENT" != "0" && "$FORCE_ALIGNMENT" != "false" && -e $ALIGNER_OUTPUT_FILE ]]; then
    echo "forcing overwrite of alignment file: "`basename $ALIGNER_OUTPUT_FILE`
    rm -f $ALIGNER_OUTPUT_FILE
fi
if [ -e $ALIGNER_OUTPUT_FILE ]; then
    echo "alignment file already exists: "`basename $ALIGNER_OUTPUT_FILE`
    echo "manually remove it or set --force-alignment to allow alignment to run again"
else

#------------------------------------------------------------------
# provide log feedback
#------------------------------------------------------------------
echo "aligning read sequences to $GENOME using $ALIGNER $FASTP_MERGE_NOTICE"
echo "  platform: $SEQUENCING_PLATFORM"
echo "  pair type: $READ_PAIR_TYPE"
echo "  length type: $READ_LENGTH_TYPE"
if [ "$IS_FIXED_READ_LENGTH" == "" ]; then
    echo "  fixed length: FALSE"
else
    echo "  fixed length: TRUE"
    echo "  read length: $FIXED_READ_LENGTH"
fi
if [ "$ADAPTER_SEQUENCE" != "" ]; then
    echo "  trim adapter r1: $ADAPTER_SEQUENCE"
fi
if [ "$ADAPTER_SEQUENCE_R2" != "" ]; then
    echo "  trim adapter r2: $ADAPTER_SEQUENCE_R2"
fi
echo "  min insert: $MIN_INSERT_SIZE"
echo "  bandwidth: $BANDWIDTH" 
echo "  aln mode: $ALIGNMENT_MODE" 
echo "  reference genome: $GENOME_FASTA"
echo "  genotype: $GENOTYPE_DIR"
echo "  read 1 files:"
echo "$READ_1_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
if [ "$READ_PAIR_TYPE" = "paired" ]; then
    echo "  read 2 files:"
    echo "$READ_2_FILES" | perl -ne 'print "    ".join("\n    ", split(" ", $_)),"\n"'
fi
echo "  output: $ALIGNER_OUTPUT_FILE"

#------------------------------------------------------------------
# platform-specific fastp pre-alignment read processing with:
#       adapter trimming
#       read merging
#       aggressive quality filtering
#   large numbers of threads do not improve fastp speed (its already fast)
#   dup evaluation is memory intensive and not needed
#   as desired, fastp prefers read1 bases over read2 bases when merging
#------------------------------------------------------------------
export FASTP_COMMAND="fastp \
    --thread 3 \
    --stdin --stdout $FASTP_INTERLEAVED_OPTIONS \
    --dont_eval_duplication \
    --length_required $MIN_INSERT_SIZE $FASTP_MERGE_OPTIONS $FASTP_ADAPTER_OPTIONS $FASTP_POLY_X_OPTIONS \
    --qualified_quality_phred $QUALIFIED_QUALITY_PHRED --unqualified_percent_limit $UNQUALIFIED_PERCENT_LIMIT \
    --average_qual $AVERAGE_QUAL \
    --html $FASTP_LOG_PREFIX.html --json $FASTP_LOG_PREFIX.json \
    --report_title \"$DATA_NAME\" 2>$FASTP_LOG_PREFIX.log.txt"
    # --cut_right --cut_window_size $CUT_WINDOW_SIZE --cut_mean_quality $CUT_MEAN_QUALITY \
export FASTP_BUFFER=${TMP_DIR_WRK_SMALL}/FASTP_BUFFER.fastq

#------------------------------------------------------------------
# reference-specific preparation (cannot export variables, do that in Workflow.sh)
#------------------------------------------------------------------
source $ACTION_DIR/prepare.sh

#------------------------------------------------------------------
# process reads and align to genome
#------------------------------------------------------------------
echo
echo "aligning read sequences"
export -f checkPipe

# pull reads from various source types to a consistent (interleaved) format
perl $PIPELINE_SHARED/prepare_fastq.pl |

# use fastp for one-pass adapter trimming, quality filtering, and read merging
#   use buffered fastp until its memory leak is fixed: https://github.com/OpenGene/fastp/issues/392
#   the leak was supposedly fixed in fastp v0.24 but memory climbing still occurs
#   continue using buffered_fastp.pl as it adds little overhead since aligner is the slow step
perl $PIPELINE_SHARED/buffered_fastp.pl | 

# tweak the way read pair merge status is reported in QNAME line
# report summary of fastp actions
perl $PIPELINE_SHARED/adjust_merge_tags.pl |

# run the alignment and any required post-processing
# handled by reference-type-specific alignment script
bash $ALIGN_SCRIPT
checkPipe

echo "done"

fi


# $ fastp --help
# usage: fastp [options] ... 
# options:
#   -i, --in1                            read1 input file name (string [=])
#   -o, --out1                           read1 output file name (string [=])
#   -I, --in2                            read2 input file name (string [=])
#   -O, --out2                           read2 output file name (string [=])
#       --unpaired1                      for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it. (string [=])
#       --unpaired2                      for PE input, if read2 passed QC but read1 not, it will be written to unpaired2. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])
#       --overlapped_out                 for each read pair, output the overlapped region if it has no any mismatched base. (string [=])
#       --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
#   -m, --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
#       --merged_out                     in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=])
#       --include_unmerged               in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.
#   -6, --phred64                        indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
#   -z, --compression                    compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
#       --stdin                          input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
#       --stdout                         stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.
#       --interleaved_in                 indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
#       --reads_to_process               specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
#       --dont_overwrite                 don't overwrite existing files. Overwritting is allowed by default.
#       --fix_mgi_id                     the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.
#   -V, --verbose                        output verbose log information (i.e. when every 1M reads are processed).
#   -A, --disable_adapter_trimming       adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
#   -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
#       --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])
#       --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
#       --detect_adapter_for_pe          by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
#   -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
#   -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
#   -b, --max_len1                       if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
#   -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
#   -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
#   -B, --max_len2                       if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
#   -D, --dedup                          enable deduplication to drop the duplicated reads/pairs
#       --dup_calc_accuracy              accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
#       --dont_eval_duplication          don't evaluate duplication rate to save time and use less memory.
#   -g, --trim_poly_g                    force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
#       --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
#   -G, --disable_trim_poly_g            disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
#   -x, --trim_poly_x                    enable polyX trimming in 3' ends.
#       --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
#   -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
#   -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
#   -r, --cut_right                      move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
#   -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
#   -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
#       --cut_front_window_size          the window size option of cut_front, default to cut_window_size if not specified (int [=4])
#       --cut_front_mean_quality         the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
#       --cut_tail_window_size           the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
#       --cut_tail_mean_quality          the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
#       --cut_right_window_size          the window size option of cut_right, default to cut_window_size if not specified (int [=4])
#       --cut_right_mean_quality         the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
#   -Q, --disable_quality_filtering      quality filtering is enabled by default. If this option is specified, quality filtering is disabled
#   -q, --qualified_quality_phred        the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
#   -u, --unqualified_percent_limit      how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
#   -n, --n_base_limit                   if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
#   -e, --average_qual                   if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
#   -L, --disable_length_filtering       length filtering is enabled by default. If this option is specified, length filtering is disabled
#   -l, --length_required                reads shorter than length_required will be discarded, default is 15. (int [=15])
#       --length_limit                   reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
#   -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
#   -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
#       --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
#       --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
#       --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
#   -c, --correction                     enable base correction in overlapped regions (only for PE data), default is disabled
#       --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
#       --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
#       --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
#   -U, --umi                            enable unique molecular identifier (UMI) preprocessing
#       --umi_loc                        specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
#       --umi_len                        if the UMI is in read1/read2, its length should be provided (int [=0])
#       --umi_prefix                     if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
#       --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])
#       --umi_delim                      delimiter to use between the read name and the UMI, default is : (string [=:])
#   -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
#   -P, --overrepresentation_sampling    one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])
#   -j, --json                           the json format report file name (string [=fastp.json])
#   -h, --html                           the html format report file name (string [=fastp.html])
#   -R, --report_title                   should be quoted with ' or ", default is "fastp report" (string [=fastp report])
#   -w, --thread                         worker thread number, default is 3 (int [=3])
#   -s, --split                          split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
#   -S, --split_by_lines                 split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
#   -d, --split_prefix_digits            the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
#       --cut_by_quality5                DEPRECATED, use --cut_front instead.
#       --cut_by_quality3                DEPRECATED, use --cut_tail instead.
#       --cut_by_quality_aggressive      DEPRECATED, use --cut_right instead.
#       --discard_unmerged               DEPRECATED, no effect now, see the introduction for merging.
#   -?, --help                           print this message


# on AvitiPE, must trim adapters by sequence (not just overlap) since short inserts will sequence all the way around the circle
# this leads fastp to merge the adapters, not the reads...

# @2581:0:0:0
# gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg sssssssssssss yyyyyyyyyyyyyyyyyyyy iiiiiiiii 
# CTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAG AGATCGGAAGAGC ACACGTCTGAACTCCAGTCA CCCGCGGTT ATCTCGTATGCCGTCTTCTGCTTGATACGGCGACCACCGAGATCTACA
#                                                                                                           ================================================
#                                                                                                                                 555555555555555555
#                                                                                                           777777777777777777777777
# +
# GLLLLLLLLMMMMMMNMNNNMNNNNNNNNNNNNNNNNNNMNNLNNNKNNNNNNNNLNNNMNNNNMNINNMNNNNNNNNNMNNNNNNNNNLNNNNNNKNMMNNLN4NINNMNLMNLLMMNMNMMNMMLJGMLMLLLMM<LMLLLMLIKGJF
# @2581:0:0:0
# CTCACTGCAACGTCTGCCTCCCAGGCTCAAGCGATTCTCCTGCCTCAGACTCCAGAGTAG AGATCGGAAGAGC GTCGTGTAGGGAAAGAGTGT CTAGCGCTG TGTAGATCTCGGTGGTCGCCGTATCAAGCAGAAGACGGCATACGAGAT
#                                                                                                           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                                                                                                                   777777777777777777777777
#                                                                                                                   555555555555555555(rc)
# +
# FJHHIIIIMMMMMMMNNNNNNNNMNNNNNNNNNNNNNNNNNNNNMNNNNNNNNNNNNNNMNNNNNNMNNMNNNNNNMNMNNMMNNNNNNMNNNMNNNNMNMMNMNMMMMMMMMJMMMMMMMMMLLMJLMMMMIMLKMLLLLKKLLKLLKI
# RC of read1: CTCACTGCAACGTCTGCCTCCCAGGCTCAAGCGATTCTCCTGCCTCAGACTCCAGAGTAG
# Illumina adapter per fastp: AGATCGGAAGAGC ACACGTCTGAACTCCAGTCA
#                             AGATCGGAAGAGC GTCGTGTAGGGAAAGAGTGT
# and its rc: TGACTGGAGTTCAGACGTGTG CTCTTCCGATCT
# TGTAGATCTCGGTGGTCGCCGTATCAAGCAGAAGACGGCATACGAGAT AACCGCGGGTGACTGGAGTTCAGACGTGT
# P5: 5' AATGATACGGCGACCACCGA 3' rc TCGGTGGTCGCCGTATCATT
# P7: 5' CAAGCAGAAGACGGCATACGAGAT 3'  rc ATCTCGTATGCCGTCTTCTGCTTG

# @2581:0:0:0 merged_150_102
# CTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTGATACGGCGACCACCGAGATCTACA
# CAGCGCTAGACACTCTTTCCCTACACGACGCTCTTCCGATCTCTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAG
# +
# GLLLLLLLLMMMMMMNMNNNMNNNNNNNNNNNNNNNNNNMNNLNNNKNNNNNNNNLNNNMNNNNMNINNMNNNNNNNNNMNNNNNNNNNLNNNNNNKNMMNNLN4NINNMNLMNLLMMNMNMMNMMLJGMLMLLLMM<LMLLLMLIKGJFMMNMNNNNMNNNMNNNNNNMMNNMNMNNNNNNMNNMNNNNNNMNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNMNNNNNNNNMMMMMMMIIIIHHJF
# CTCACTGCAACGTCTGCCTCCCAGGCTCAAGCGATTCTCCTGCCTCAGACTCCAGAGTAG AGATCGGAAGAGC GTCGTGTAGGGAAAGAGTGTCTAGCGCTG

# @2581:0:0:0:2:150:102
# CTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTGATACGGCGACCACCGAGATCTACACAGCGCTAGACACTCTTTCCCTACACGACGCTCTTCCGATCTCTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAG
# +
# GLLLLLLLLMMMMMMNMNNNMNNNNNNNNNNNNNNNNNNMNNLNNNKNNNNNNNNLNNNMNNNNMNINNMNNNNNNNNNMNNNNNNNNNLNNNNNNKNMMNNLN4NINNMNLMNLLMMNMNMMNMMLJGMLMLLLMM<LMLLLMLIKGJFMMNMNNNNMNNNMNNNNNNMMNNMNMNNNNNNMNNMNNNNNNMNNNNNNNNNNNNNNMNNNNNNNNNNNNNNNNNNNNMNNNNNNNNMMMMMMMIIIIHHJF

# the problem resolves if adapter sequences are provided
# @2581:0:0:0 merged_60_0
# CTACTCTGGAGTCTGAGGCAGGAGAATCGCTTGAGCCTGGGAGGCAGACGTTGCAGTGAG
# +
# GLLLLLLLLMMMMMMNMNNNMNNNNNNNNNNNNNNNNNNMNNLNNNKNNNNNNNNLNNNM
