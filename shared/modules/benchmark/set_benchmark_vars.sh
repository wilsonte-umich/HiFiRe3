# action:
#     set environment variables for performing benchmarking
# requires:
#     source $MODULES_DIR/library/set_library_vars.sh
# usage:
#     source $MODULES_DIR/analyze/set_analysis_vars.sh

# analysis identifiers
export BENCHMARK_IDENTIFIER=${DATA_NAME}.benchmark
export BENCHMARK_PREFIX=${TASK_DIR}/${BENCHMARK_IDENTIFIER}

# SV benchmark output
export BENCHMARK_BAM=${BENCHMARK_PREFIX}.coordinate.bam
export SNIFFLES_VCF=${BENCHMARK_PREFIX}.sniffles.vcf.gz
export SNIFFLES_FRAGSIZES_VCF=${BENCHMARK_PREFIX}.sniffles.frag_sizes.vcf.gz
export ADMIXED_CALLS_CSV=${BENCHMARK_PREFIX}.admixed_calls.csv
