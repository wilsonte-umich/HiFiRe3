use strict;
use warnings;

# action:
#   check whether alignments pass alignment-level quality filters:
#       all alignments, even single, non-SV alignments per read:
#           minimum MAPQ
#           maximium gap-corrected base divergence
#       alignment in SV reads only:
#           minimum alignment length, measured as reference base span
#           (ONT only) minimum average base quality across alignment
#               with Dorado, individual base qualities range from Q0 to Q50
#               see script output for the range of average base qualities observed over SV flanking alignments
#   when an alignment fails, remove it and all alignments further 3' on the read
#       this action suppresses untrusted SV junctions that would otherwise include the failed alignment
# input:
#   name-sorted SITE_SAM on STDIN, ordered from read 5' to 3' end
# output: 
#   name-sorted SITE_SAM on STDOUT
#   where this script removes some alignments and renumbers some reads

# initialize reporting
our $script = "enforce_quality_filters";
our $error  = "$script error";
my ($nSequences, $nReadsIn, $nSvReadsIn, $nReadsOut, $nSvReadsOut) = (0) x 10;
my (@nAlnsRejected, %nAlnsRejected, @avgBaseQual);

# load dependencies
my $perlUtilDir = "$ENV{MODULES_DIR}/utilities/perl";
map { require "$perlUtilDir/$_.pl" } qw(workflow numeric);
map { require "$perlUtilDir/sequence/$_.pl" } qw(general);
resetCountFile();

# environment variables
fillEnvVar(\our $MIN_MAPQ,          'MIN_MAPQ');
fillEnvVar(\our $MAX_DIVERGENCE,    'MAX_DIVERGENCE');
fillEnvVar(\our $MIN_FLANK_LEN,     'MIN_FLANK_LEN');
fillEnvVar(\our $MIN_AVG_BASE_QUAL, 'MIN_AVG_BASE_QUAL');
fillEnvVar(\our $HAS_BASE_ACCURACY, 'HAS_BASE_ACCURACY');

# constants
use constant {
    S_QNAME             => 0, # SITE_SAM fields
    S_FLAG              => 1,
    S_RNAME             => 2,
    S_POS1              => 3, # 1-based
    S_MAPQ              => 4,
    S_CIGAR             => 5,
    CH_TAG              => 6,
    TL_TAG              => 7,
    DE_TAG              => 8,
    HV_TAG              => 9,
    N_REF_BASES         => 10,
    N_READ_BASES        => 11,
    BLOCK_N             => 12,
    SITE_INDEX1_1       => 13,
    SITE_POS1_1         => 14,
    SITE_DIST_1         => 15,
    SITE_INDEX1_2       => 16,
    SITE_POS1_2         => 17,
    SITE_DIST_2         => 18,
    SEQ_SITE_INDEX1_2   => 19,
    SEQ_SITE_POS1_2     => 20,
    IS_END_TO_END       => 21,
    READ_HAS_JXN        => 22,
    TARGET_CLASS        => 23,
    S_SEQ               => 24,
    S_QUAL              => 25,
    CS_TAG              => 26,
    #-------------
    SPLIT_TO_S_QUAL     => 27,
    #-------------
    _IS_PAIRED      => 1, # SAM FLAG bits
    _PROPER_PAIR    => 2,
    _UNMAPPED       => 4,
    _MATE_UNMAPPED  => 8,
    _REVERSE        => 16,
    _MATE_REVERSE   => 32,
    _FIRST_IN_PAIR  => 64,
    _SECOND_IN_PAIR => 128,
    _SECONDARY      => 256,
    _FAILED_QC      => 512,
    _DUPLICATE      => 1024,
    _SUPPLEMENTAL   => 2048,
    # -------------
    REJECT_MAPQ             => "mapq", # these are just hash keys, not report strings
    REJECT_DIVERGENCE       => "divergence",
    REJECT_FLANK_LEN        => "flank",
    REJECT_AVG_BASE_QUAL    => "baseQual",
    #-------------
    AVG_BASE_QUAL_BIN_SIZE => 5,
};

# parse input SAM
my ($prevQName, @alns);
while(my $aln = <STDIN>){
    chomp $aln;
    my @aln = split("\t", $aln, SPLIT_TO_S_QUAL);
    if($prevQName and $prevQName ne $aln[S_QNAME]){
        processQName();
        @alns = ();
    }
    push @alns, \@aln;
    $prevQName = $aln[S_QNAME];
}
processQName(); 

# print summary information
printCount(commify($nSequences),    'nSequences',   'input sequences');
printCount(commify($nReadsIn),      'nReadsIn',     'input reads');
printCount(commify($nSvReadsIn),    'nSvReadsIn',   'input SV reads');
printCount(commify($nReadsOut),     'nReadsOut',    'output reads');
printCount(commify($nSvReadsOut),   'nSvReadsOut',  'output SV reads');
printCount(commify($nAlnsRejected{REJECT_MAPQ}          || 0), 'nRejectMapq',           "reads were truncated with alignment MAPQ < $MIN_MAPQ");
printCount(commify($nAlnsRejected{REJECT_DIVERGENCE}    || 0), 'nRejectDivergence',     "reads were truncated with alignment gap-corrected divergence > $MAX_DIVERGENCE");
printCount(commify($nAlnsRejected{REJECT_FLANK_LEN}     || 0), 'nRejectFlankLen',       "SV reads were truncated with alignment length < $MIN_FLANK_LEN");
printCount(commify($nAlnsRejected{REJECT_AVG_BASE_QUAL} || 0), 'nRejectAvgBaseQual',    "SV reads were truncated with alignment avg. base QUAL < $MIN_AVG_BASE_QUAL");
foreach my $i(0..$#nAlnsRejected){
    printCount(commify($nAlnsRejected[$i] || 0), 'nReject_'.($i + 1), "reads were truncated at alignment".($i + 1));
}
if(!$HAS_BASE_ACCURACY){
    print STDERR "SV alignment avgBaseQuals\n";
    for (my $avgBaseQual = 0; $avgBaseQual <= 50; $avgBaseQual += AVG_BASE_QUAL_BIN_SIZE){
        my $bin = $avgBaseQual / AVG_BASE_QUAL_BIN_SIZE;
        print STDERR join("\t", $avgBaseQual, $avgBaseQual[$bin] || 0), "\n";
    }
}

# process QNAME alignment groups
sub processQName {
    $nSequences++;

    # check alignment-level quality metrics
    # reject everything distal to a rejected alignment, inclusive
    my @outAlns;
    my $readHasSv = @alns > 1;
    $readHasSv and $nSvReadsIn++;
    foreach my $i(0..$#alns){
        my $aln = $alns[$i];

        # rejection criteria are enforced sequentially in order of efficiency
        # i.e., frequent, easy rejections are checked first

        # criteria enforced on all alignments, even single, non-SV alignments
        if($$aln[S_MAPQ] < $MIN_MAPQ){
            $nAlnsRejected{REJECT_MAPQ}++;
            $nAlnsRejected[$i]++;
            last;
        }
        if($$aln[DE_TAG] > $MAX_DIVERGENCE){
            $nAlnsRejected{REJECT_DIVERGENCE}++;
            $nAlnsRejected[$i]++;
            last;
        }

        # criteria only enforced when reads have SV junctions, i.e., multiple alignments
        if($readHasSv){
            if(getEnd(@$aln[S_POS1, S_CIGAR]) - $$aln[S_POS1] < $MIN_FLANK_LEN){
                $nAlnsRejected{REJECT_FLANK_LEN}++;
                $nAlnsRejected[$i]++;
                last;
            }
            if(
                !$HAS_BASE_ACCURACY and # thus, this slow check only performed on ONT or other low accuracy platform
                $$aln[S_QUAL] ne "*"
            ){
                my $avgBaseQual = getAvgQual($$aln[S_QUAL]);
                $avgBaseQual[int($avgBaseQual / AVG_BASE_QUAL_BIN_SIZE + 0.5)]++;
                if($avgBaseQual < $MIN_AVG_BASE_QUAL){
                    $nAlnsRejected{REJECT_AVG_BASE_QUAL}++;
                    $nAlnsRejected[$i]++;
                    last;
                }
            }
        }
        push @outAlns, $aln;
    }

    # print the final kept alignments, if any
    @outAlns and $nReadsOut++;
    @outAlns > 1 and $nSvReadsOut++;
    foreach my $aln(@outAlns){
        print join("\t", @$aln), "\n";
    }
}

1;

# example of a ubam QUAL after running Dorado via `basecall ONT`, but before minimap2
# maximum base quality appears to be S, i.e., Phred Q50
# very low values are also present, e.g., # == Q2, % = Q4

#;:;==FEGGHJB@6669=A=0HSSHJIHGGGJGMPSSSSKSKHLJLSP>KSLNGGKHHJHSSSSSSSSSNNSSGKHSSNKJLSGHHGFFFIG@?@@@@<H8877;:DFGEKJHIIHLSSSSJFFF?=:99LIHIFCBEISSSSSSSSSSSJFEB43333LSNQSSLJSKMILSSSSSSSSSSKSSSSSSSMSSSSRSFFSSQLJIHHJ444SCCDKPNSSSLSSSSPJM==<==LJQSSSSSSSKSSSSPSSSSSMEFHMFJD654211111211>=??@?EJRLDFHHSDDSSSHJLDFFSGSSRIJISSKSLSHHFFAB@@@@SSSSMJRKQSGKHDBBDDFGJKIIHPISSSSHGSJG:::::SOSLHLSEEFFGA@A@ELHJSCOIKEGCBB?>>??SNLSSIMKKSSSJJLIGIKGD9865588KIGGSSSSSSSGHGHILSSGGFEFGFGGGSSNMSHIKQQJGGFDDFFGIJOOOSHGDECHJLHIHDA??>EFGHCBCE887700008007721(((((-.:>>>2222CCPIIIGGJLJLRIIMSOC))))CIHHKPEKNHGIGHINMM999FFSEFF6==211137222334445:@A@?56@GGHILKIOFCHLA66DDKSMSSOKSIGFSSSSOLSPMLEDIGGSSSSSSSSNMQ?????SSOPKSLMKMNFSSSLSSOSQLHSSNLSPNSSSSSSSKSSSSKSNJSSGSDBC9,*.,.)(()11(((22<11>CCGIISSSSNLLSKIJFGHGC@@66HGGGSSLSROIKNIE@@DG@@@@@@HHGGEEDHKFFFHGB><887B@DFFGFSNFFHFFISJHKEGHHEDDDEHJJOIFFFFIIIKGGGGSSSSSHNKISMKMGGGHFBCBACJSIFISKSLFFFFFDGEEDCGCFGCDDDCACBBCFFFDCC@?=:3>AEBBDBCCEIGKOLMSLHBFSRMKIEGDEGIHSSPSSC=?>CEFIGBBAABHQSMIJNHIBB;::::88888BIJQL<;;;<SNMJLH=><<<AAGMNPMA?BA<4JSOLKLHLNQSISKHJFAB=:5//...-.-4/,)'&%()13<?A6BCBDDCCCAAAACGHFFCCA@A@@@B@?@@ADC7532357?FDECABAAAA@>>??@CBBBCCCCCCCDDDCCDFFGGFCCCCCCCCCC>;;;:99346622...//////GECCDCDDDCFDCCAAACAAABACCB@@A@@@>>?;:90---------<BA@ACD>=>=<=?ABCCC>??>>@?==??@@@???@@@A@=9999:CA===CACBCCABB::?=1+==98976553333;A@:::BCGCACBCCCGHGFEE@?>==?@ACD??>>=>>>@

# example avgBaseQual on ~20K alignment from ONT SV reads
# avgBaseQual	Count
# 20	0
# 21	0
# 22	0
# 23	0
# 24	1
# 25	1
# 26	8
# 27	14
# 28	38
# 29	92
# 30	286
# 31	495
# 32	1010
# 33	1181
# 34	1005
#------------
# 35	1502
# 36	1329
# 37	1080
# 38	925
# 39	825
# 40	993
# 41	1018
# 42	1479
# 43	2084
# 44	2590
# 45	1677
# 46	312
# 47	10
# 48	0
# 49	0
# 50	0
