#!/bin/bash

# common UCSC URLs
UCSC_GOLDEN_PATH=https://hgdownload.soe.ucsc.edu/goldenPath
HG38_ANALYSIS_SET="${UCSC_GOLDEN_PATH}/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"

# shared download action support
download_and_unzip_file () {
    if [ ! -f $TO_FILE ]; then
        if [ ! -f $GZ_FILE ]; then
            echo "  downloading $FILE_DESCRIPTION"
            echo "    from: $URL"
            echo "    to:   $GZ_FILE"
            wget -O $GZ_FILE $URL
        fi
        if [ "$LEAVE_ZIPPED" == "" ]; then
            echo "  unzipping $FILE_DESCRIPTION"
            gunzip $GZ_FILE
        fi
    else 
        echo "  already exists: $FILE_DESCRIPTION"
    fi    
}
download_and_unzip_UCSC () {
    UCSC_FILE=$1
    TO_FILE=$2
    FILE_DESCRIPTION=$3
    LEAVE_ZIPPED=$4
    if [ "$LEAVE_ZIPPED" == "" ]; then
        GZ_FILE=${TO_FILE}.gz
    else 
        GZ_FILE=${TO_FILE}
    fi
    URL=${UCSC_GOLDEN_PATH}/${GENOME}/${UCSC_FILE}
    if [ "$FORCE_URL" != "" ]; then
        URL=${FORCE_URL}
        FORCE_URL=""
    fi
    download_and_unzip_file
}
download_and_index_genome () {
    FORCE_URL=$1 # to support downloading an analysis set instead of just the genome
    download_and_unzip_UCSC \
        bigZips/${GENOME}.fa.gz \
        ${GENOME_FASTA} \
        "genome fasta"
    if [ ! -f $GENOME_FASTA.fai ]; then
        if [ "$APPEND_EBV" != "" ]; then
            FORCE_URL=${HG38_ANALYSIS_SET}
            WRK_DIR=`dirname ${GENOME_FASTA}`
            HG38_FILE=${WRK_DIR}/hg38.fa
            EBV_FILE=${WRK_DIR}/chrEBV.fa
            if [ ! -f ${EBV_FILE} ]; then # get the EBV genome from the hg38 analysis set
                download_and_unzip_UCSC \
                    NA \
                    ${HG38_FILE} \
                    "hg38 analysis set"
                echo "  indexing analysis set"
                samtools faidx ${HG38_FILE}
                echo "  extracting chrEBV"
                samtools faidx ${HG38_FILE} chrEBV > ${EBV_FILE}
                rm ${HG38_FILE}
                rm ${HG38_FILE}.fai # leave chrEBV.fa but remove the larger, now obsolete hg38 file
            fi
            echo "  appending chrEBV"
            cat ${EBV_FILE} >> ${GENOME_FASTA}
        fi
        echo "  indexing genome fasta"
        samtools faidx $GENOME_FASTA
    else 
        echo "  already exists: genome fasta index"
    fi
}
download_genome_gaps () {
    download_and_unzip_UCSC \
        database/gap.txt.gz \
        ${GENOME_GAPS_FILE} \
        "gaps file"
}
download_ENCODE_exclusions () {
    TO_FILE=${GENOME_EXCLUSIONS_BED}
    GZ_FILE=${TO_FILE}.gz
    FILE_DESCRIPTION="excluded regions"
    LEAVE_ZIPPED=""
    URL=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/${GENOME}-blacklist.v2.bed.gz
    download_and_unzip_file
}
download_excluderegions () { # for newer genomes lacking Boyle lab exclusions
    REMOTE_FILE=$1
    TO_FILE=${GENOME_EXCLUSIONS_BED}
    if [ ! -f $TO_FILE ]; then
        echo
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "at present, the best choice for obtaining exclusion regions for genome ${GENOME}"
        echo "is to manually download the bed format (not bed.gz format) file:"
        echo "    $REMOTE_FILE"
        echo "from Google Drive:"
        echo "    https://drive.google.com/drive/folders/1sF9m8Y3eZouTZ3IEEywjs2kfHOWFBSJT/"
        echo "to this exact path and file name:"
        echo "    $TO_FILE"
        echo "for details on the source of these files, see:"
        echo "    https://github.com/dozmorovlab/excluderanges"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
    fi
}
download_annotation_gtf () {
    TO_FILE=${ANNOTATION_GTF}
    GZ_FILE=${TO_FILE}
    FILE_DESCRIPTION="genome GTF"
    LEAVE_ZIPPED="TRUE"
    URL=${UCSC_GOLDEN_PATH}/${GENOME}/bigZips/genes/${GENOME}.ncbiRefSeq.gtf.gz
    download_and_unzip_file

    if [ ! -f $GENES_BED ]; then
        echo "  converting GTF to genes BED"
        zcat ${ANNOTATION_GTF} |
        awk '$3 == "exon" && $1 !~ /_/' | # skip alt contigs, which leads to gene name redundancy
        perl -ne '
            chomp;
            my @f = split /\t/;
            $f[8] =~ m/gene_name "(.+?)"/ and print join("\t",
                $f[0],     # seqname
                $1,        # gene_name
                $f[3] - 1, # start0
                $f[4],     # end1
                $f[6]      # strand
            ), "\n";
        ' | 
        sort -k1,2 | # allow same gene names on two chromosomes; required for chrX/Y pseudoautosomal regions
        bedtools groupby -g 1,2 -c 3,4,5 -o min,max,mode  | # thus, the span of all exons in all transcripts of a named gene
        awk '
            BEGIN { OFS="\t" }
            { print $1, $3, $4, $2, $4 - $3, $5 }
        ' |
        awk '$5 < 5000000' | # skip rare gene names dispersed too widely on a chrom
        sort -k1,1 -k2,2n -k3,3n |
        gzip -c > ${GENES_BED} # output has mostly one bed row per gene, but gene names are _not_ universally unique, see above
    fi
}

# execute download actions customized for each specific supported genome
echo
echo "downloading resource files for reference genome ${GENOME}"
if [ "$GENOME" == "hs1" ]; then
    APPEND_EBV="TRUE"
    download_and_index_genome ""
    echo "  touching genome gaps (hs1/CHM13 has no gaps!)"
    touch ${GENOME_GAPS_FILE} # CHM13 has no gaps! create the empty file anyway for consistency
    download_excluderegions T2T.excluderanges.bed
    download_annotation_gtf

elif [ "$GENOME" == "hg38" ]; then
    APPEND_EBV=""
    download_and_index_genome ${HG38_ANALYSIS_SET}
    download_genome_gaps
    download_ENCODE_exclusions
    download_annotation_gtf

elif [ "$GENOME" == "mm39" ]; then
    APPEND_EBV=""
    download_and_index_genome ""
    download_genome_gaps
    download_excluderegions mm39.excluderanges.bed
    download_annotation_gtf

elif [ "$GENOME" == "dm6" ]; then
    APPEND_EBV=""
    download_and_index_genome ""
    download_genome_gaps
    download_ENCODE_exclusions
    download_annotation_gtf
fi

echo
echo "done"
