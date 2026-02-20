#!/bin/bash

# for all files, concatenate files in order GENOME1 then GENOME2
# appending _$GENOME1 or _$GENOME2 to the chromosome name column

# create the app composite genome metadata file
echo "creating app metadata file"
perl ${ACTION_DIR}/create_app_metadata.pl

# concatenate genome gaps
if [ ! -f $GENOME_GAPS_FILE ]; then
    echo "concatenating genome gaps"
    cat \
        <(awk 'BEGIN{OFS="\t"}{$2 = $2"_'$GENOME1'"; print $0 }' $GENOME1_GAPS_FILE) \
        <(awk 'BEGIN{OFS="\t"}{$2 = $2"_'$GENOME2'"; print $0 }' $GENOME2_GAPS_FILE) \
        > $GENOME_GAPS_FILE
fi

# concatenate genome exclusions
if [ ! -f $GENOME_EXCLUSIONS_BED ]; then
    echo "concatenating exclusions bed"
    cat \
        <(awk 'BEGIN{OFS="\t"}{$1 = $1"_'$GENOME1'"; print $0 }' $GENOME1_EXCLUSIONS_BED) \
        <(awk 'BEGIN{OFS="\t"}{$1 = $1"_'$GENOME2'"; print $0 }' $GENOME2_EXCLUSIONS_BED) \
        > $GENOME_EXCLUSIONS_BED
fi

# concatenate gene annotations
if [ ! -f $GENES_BED ]; then
    echo "concatenating genes bed"
    cat \
        <(zcat $GENES1_BED | awk 'BEGIN{OFS="\t"}{$1 = $1"_'$GENOME1'"; print $0 }') \
        <(zcat $GENES2_BED | awk 'BEGIN{OFS="\t"}{$1 = $1"_'$GENOME2'"; print $0 }') |
        gzip -c > $GENES_BED
fi

# concatenate genome fasta
if [ ! -f $GENOME_FASTA ]; then
    echo "concatenating genome fasta"
    cat \
        <(awk 'BEGIN{OFS="\t"}{if($1 ~ /^>/) print $1"_'$GENOME1'"; else print $0 }' $GENOME1_FASTA) \
        <(awk 'BEGIN{OFS="\t"}{if($1 ~ /^>/) print $1"_'$GENOME2'"; else print $0 }' $GENOME2_FASTA) \
        > $GENOME_FASTA

    echo "indexing genome fasta"
    samtools faidx $GENOME_FASTA
fi
