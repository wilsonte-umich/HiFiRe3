
echo "determining RE frag sizes from read consensuses"

${HF3_TOOLS_BIN} digest_sv_alleles > ${SNIFFLES_FRAGSIZES_VCF}
checkPipe

echo "done"
