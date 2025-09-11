use strict;
use warnings;

# action:
#   describe the match/SNV/indel combinations of the 5'-most alignment of all reads
#   use an efficient fragment index to count and sort unique combinations
#   work is parallelized by chromosome for speed gains and sorting efficiency
#   this script is only applicable to sequencing platforms where HAS_BASE_ACCURACY=TRUE
#   match/SNV/indel combinations in an alignment are described by:
#       readOffset5 = number of bases to place the alignment onto the read relative to its 5' end
#                     thus, the number of bases clipped from the read 5' end upstream of the first vaTag : operation
#       vaTag = reformatted cs:Z: tag for the 5'-most alignment of a read, in reference strand orientation
#           see va_tag_utils.pl for a detailed description of the va tag
#           all vaTag coordinate and base _values_ are relative to the top genome strand, but _ordered_ relative to the alignment strand
#       qsTag = tag that carries boolean passed-quality operations in parallel to vaTag, where qsTag value are:
#               vaTag : operations = always 1 = passed quality (must be high or base alignment would have failed)
#               vaTag alignment length  > 0 = 0|1 indicating if the average quality of the variant read bases is >= MIN_SNV_INDEL_QUAL
#               vaTag alignment length == 0 = 0|1 indicating if the average quality of the two bases flanking the isolated deletion is >= MIN_SNV_INDEL_QUAL
#           qsTag operation values are overridden to 1 if readHasSnv is FALSE, i.e., all read SNVs or indels matched a known haplotype
#       xRead1 = for merged PE, number of read 5' bases derived from read1; for others, 0|1 indicating if this was read1
#       xRead2 = for merged PE, number of read 3' bases derived from read2; for others, 0|1 indicating if this was read2
#           for SE or orphan PE, xRead1 is always 1, xRead2 is always 0
#           for unmerged PE, the index does not yet aggregate read1 and read2 (since read1 typically has higher accuracy)
#   to facilitate construction of base pileups downstream, each reference strand is indexed separately:
#        in the top-strand output file, strandIndex0=0:
#            alignments are sorted in genome order
#            va tag operations are in 5'-3' order on reference
#            base coordinates used for sorting are the same as the reference genome
#            leftPos1 is the 5'-most base of the read alignment, leftmost on reference
#        in the bottom-strand output file, strandIndex0=1:
#            alignments are sorted in reverse genome order
#            va tag operations are reversed to be in 3'-5' order on reference (5'-3' on reversed reference)
#            base coordinates used for sorting are chromSize - pos1 + 1 (thus, numbered 5' to 3' on the matching strand)
#            leftPos1 is also the 5'-most base of the read alignment, now rightmost on reference (leftmost on reversed reference)
#        thus, both files represent bases:
#            in 5' to 3' order on read, starting at a RE site in the case of ligFree, and 
#            in 5' to 3' order on the corresponding reference _strand_ (not always the top strand as in many other encodings)
#   a hash is used to track the unique alignments encountered for each bin as:
#       zeroPaddedLeftPos1,zeroPaddedRightPos1,readOffset5,vaTag,qsTag,xRead1,xRead2
#   where the following hash components are used for sorting but omitted from the final output:
#       leftPos1  = the lowest  reference strand coordinate of the alignment, corresponding to the first vaTag : operation
#       rightPos1 = the highest reference strand coordinate of the alignment
#           zeroPadded postions support simple and efficient lexical sorting of positions within bins
#   the indexing strategy exploits threads, loops, and numeric indices to pre-sort by:
#       chromIndex1 => strandIndex0 => binIndex0
#       where binIndex0 bins alignment 5' ends by leftPos1 for a first round of sorting into groups
#   only the 5'-most aligments of reads (not distal SV alignments) are assessed for SNVs and indels for streamlined indexing
#   indexing occurs with respect to read strand, so reads are NOT aggregated by strand
#       i.e., a variant sequenced in two strand orientations will be found on at least two indexed lines
#       no such strand aggregation is performed downstream since:
#           overlapping paired reads were previously merged to prevent redundancy in the overlap region
#           two single-end reads rarely arise from two strands of the same source duplex since PCR-free
#       thus, the potential for non-independent variant calls on two strands is considered very small
#       however, multiple indexed lines may exist for a clonal variant if sequencing errors give distinct hash keys
#   indexing is an aggregating action, so read/event identifiers are lost in the index
#   unlike for SVs, no attempt is made to retain a list to allow SNV/indels to be correlated to source reads because:
#       the main utility would be for purging duplex redundancy, and ONT is never used for SNV tracking
#       alignment-level quality filters have already been enforced
#       other essential alignment metadata are present in the index
# input:
#   SITE_SAM files generated by `analyze fragments`, one file per chromosome
# output:
#   alignment_maps = one row per unique 5'-most read alignment, two stranded files per chromosome
#   base_qualities = counts of SNP/indel base qualities stratified by SNV/indel type, useful for setting MIN_SNV_INDEL_QUAL

# initialize reporting
our $script = "index_fragments";
our $error  = "$script error";

# working variables
use vars qw($chromIndex1 $strandIndex0 @alns $varKey $snvAlnsH);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow);

# environment variables
fillEnvVar(\our $HAS_BASE_ACCURACY,     'HAS_BASE_ACCURACY');
fillEnvVar(\our $MIN_SNV_INDEL_QUAL,    'MIN_SNV_INDEL_QUAL');
fillEnvVar(\our $SNV_ALNS_PREFIX,       'SNV_ALNS_PREFIX');

# check the platform
!$HAS_BASE_ACCURACY and die "SNV/indel analysis only applies to sequencing platforms with high base accuracy\n";

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    DE_TAG              => 6,
    CS_TAG              => 7,
    XF_TAG              => 8,
    XH_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_HAPS_1         => 15,
    SITE_DIST_1         => 16,
    SITE_INDEX1_2       => 17,
    SITE_POS1_2         => 18,
    SITE_HAPS_2         => 19,
    SITE_DIST_2         => 20,
    SEQ_SITE_INDEX1_2   => 21,
    SEQ_SITE_POS1_2     => 22,
    SEQ_SITE_HAPS_2     => 23,
    IS_END_TO_END       => 24,
    READ_HAS_JXN        => 25,
    TARGET_CLASS        => 26,
    S_SEQ               => 27,
    S_QUAL              => 28,
    #-------------
    S_VA_TAG            => 29,
    S_LEFT_CLIP         => 30,
    S_RIGHT_CLIP        => 31,
    S_END_POS1          => 32,
    #-------------
    channel             => 0, # incoming qName extensions
    trim5               => 1,
    trim3               => 2,
    isMerged            => 3, # true=2 for legacy reasons
    nRead1              => 4,
    nRead2              => 5,
    splitGapReadN       => 6,
    readN               => 7,
    N_QNAME_EXTENSIONS  => 8,
    #-------------
    FALSE  => 0,
    TRUE   => 1,
    # -------------
    EVENT   => 0,
    READ1   => 1,
    READ2   => 2,
};

# create a hash key for a read and count it in the index
sub indexRead {
    my @qName = split(":", $alns[0][S_QNAME]);
    my @extensions = splice(@qName, -N_QNAME_EXTENSIONS);
    join("\t",

        # fields used for sorting, dropped from the final indexed output (the va tag carries coordinates also)
        # these positions do not include clips, they are the 5'-most and 3'-most positions of the aligned portion of the read
        sprintf("%09d", $alns[0][S_POS1]).sprintf("%09d", $alns[0][S_END_POS1]), # leftPos1.rightPos1

        # fields found in the final indexed output
        $alns[0][S_VA_TAG], # vaTag
        getQsTag($alns[0]), # qsTag
        $alns[0][S_LEFT_CLIP], # readOffset5; with S_POS1 provides an idenitifier for the 5' most end of the read
        $extensions[isMerged] ? # xRead1; these are read-level values that do not change with clipping or alignment trimming
            $extensions[nRead1] : 
            ($extensions[readN] == READ1 ? TRUE : FALSE),
        $extensions[isMerged] ? # xRead2
            $extensions[nRead2] : 
            ($extensions[readN] == READ2 ? TRUE : FALSE)
    )
}

# reverse the process to convert counted hash keys to output formats
my $varKey_start = qr/^(\d{9})\d{9}\t(.+)/;
sub printIndexedVariant {
    my ($nObserved) = @_;
    $varKey =~ m/$varKey_start/; 
    print $snvAlnsH join("\t", # file name carries chromosome and strand information
        $1 + 0, # keep alignment leftPos1 for pileup registration
        $2,
        $nObserved
    ), "\n";
}

# do the work
require "$ENV{ACTION_DIR}/index/index_frags_runner.pl";
