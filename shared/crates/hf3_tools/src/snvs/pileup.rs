use std::usize;

/// Support for making a reference pileup of error-corrected read bases.

// dependencies
use serde::Serialize;
use mdi::OutputCsv;
use super::{SnvChromWorker, SnvAnalysisTool};

/// PileupMetadata reports summary results of pileup construction.
pub struct PileupMetadata {
    pub n_chunks:          usize,
    pub n_reported_chunks: usize,
    pub n_reported_bases:  usize,
    pub reported_coverage: usize,
    pub n_variant_bases:   usize,
}
impl PileupMetadata {
    fn new() -> Self {
        PileupMetadata {
            n_chunks:          0,
            n_reported_chunks: 0,
            n_reported_bases:  0,
            reported_coverage: 0,
            n_variant_bases:   0,
        }
    }
}

/// A PileupRecord is a run of bases with the same coverage.
/// 
/// Serializable for writing to a BED format file.
#[derive(Serialize)]
struct PileupRecord<'a> {
    chrom_index: u8,
    start0:      u32, // BED half-open coordinates
    end1:        u32,
    counts:      &'a PileupCounts, // the observed counts for all bases in the run
}

/// PileupCounts holds read content counts for a single known position 
/// on a known chromosome.
#[derive(PartialEq, Eq, Clone, Serialize)]
pub struct PileupCounts { // 9 fields x 1 bytes = 9 bytes per position, so 9 * 250 Mb = 2.25 GB for human chr1
    n_ref:     u8, // to help construct identity runs, we don't store the ref base, just its count
    n_alt_a:   u8, // for memory efficiency, coverages saturate at 255
    n_alt_c:   u8,
    n_alt_g:   u8,
    n_alt_t:   u8,
    n_alt_n:   u8, // ACGT substitutions are allowed, otherwise they would be Ns
    n_alt_del: u8, // this single position is deleted in the read; implicitly allowed
    n_alt_ins: u8, // there is a read insertion following this reference base; PileupCounts does not hold the insertion
    n_alt_del_allowed: u8, // same as n_alt_del and n_alt_ins, but filtered for allowed indels only
    n_alt_ins_allowed: u8,
}
impl PileupCounts {
    /// Crete a new PileupCounts instance with all counts initialized to zero.
    fn new() -> Self {
        PileupCounts {
            n_ref:     0,
            n_alt_a:   0,
            n_alt_c:   0,
            n_alt_g:   0,
            n_alt_t:   0,
            n_alt_n:   0,
            n_alt_del: 0,
            n_alt_ins: 0,
            n_alt_del_allowed: 0,
            n_alt_ins_allowed: 0,
        }
    }

    /// Get the number of reads overlapping a genome position.
    pub fn coverage(&self) -> usize {
        self.n_ref     as usize + 
        self.n_alt_a   as usize + 
        self.n_alt_c   as usize + 
        self.n_alt_g   as usize + 
        self.n_alt_t   as usize + 
        self.n_alt_n   as usize + 
        self.n_alt_del as usize // deletion coverage reflects the number of reads that cover the variant position
    }                           // insertions are not counted as they are recorded on bases included above
}

/// ChromPileup holds a complete vector of PileupCounts for a single 
/// known chromosome at all positions.
pub struct ChromPileup(pub Vec<PileupCounts>);
impl ChromPileup {

    /// Create a new ChromPileup instance with a vector of PileupCounts 
    /// initialized to the specified capacity, i.e., chrom_size.
    pub fn with_capacity(capacity: usize) -> Self {
        ChromPileup(vec![PileupCounts::new(); capacity])
    }

    /// Increment a set of contiguous chromosome position counts for
    /// reference base(s) observed in a read alignment starting at ref_pos0.
    pub fn increment_ref(&mut self, ref_pos0: usize, len: usize) {
        for i in 0..len {
            ChromPileup::increment(&mut self.0[ref_pos0 + i].n_ref);
        }
    }

    /// Increment a chromosome position count for a specific observed 
    /// base in a read alignment at ref_pos0.
    pub fn increment_base(&mut self, ref_pos0: usize, base: char) {
        match base {
            'A' => ChromPileup::increment(&mut self.0[ref_pos0].n_alt_a),
            'C' => ChromPileup::increment(&mut self.0[ref_pos0].n_alt_c),
            'G' => ChromPileup::increment(&mut self.0[ref_pos0].n_alt_g),
            'T' => ChromPileup::increment(&mut self.0[ref_pos0].n_alt_t),
            'N' => ChromPileup::increment(&mut self.0[ref_pos0].n_alt_n),
            _  => panic!("Unexpected base character in pileup increment: {}", base),
        }
    }

    /// Increment a set of contiguous chromosome position counts for
    /// deleted base(s) observed in a read alignment starting at ref_pos0.
    pub fn increment_del(&mut self, ref_pos0: usize, len: usize, allowed: bool) {
        for i in 0..len {
            ChromPileup::increment(&mut self.0[ref_pos0 + i].n_alt_del);
            if allowed {
                ChromPileup::increment(&mut self.0[ref_pos0 + i].n_alt_del_allowed);
            }
        }
    }

    /// Increment a chromosome position count for an insertion observed 
    /// in a read alignment after ref_pos0.
    pub fn increment_ins(&mut self, ref_pos0: usize, allowed: bool) {
        ChromPileup::increment(&mut self.0[ref_pos0].n_alt_ins);
        if allowed {
            ChromPileup::increment(&mut self.0[ref_pos0].n_alt_ins_allowed);
        }
    }

    /// Increment a position count with saturation.
    fn increment(count: &mut u8) {
        *count = count.saturating_add(1);
    }

    /// Chunk a ChromPileup by runs of identical coverage counts
    /// and write to a temporary file.
    pub fn write_chunked(
        tool:   &SnvAnalysisTool,
        worker: &mut SnvChromWorker,
    ) -> PileupMetadata {
        let mut csv = OutputCsv::open_csv(
            &worker.pileup_file_path, 
            b'\t', 
            false, 
            Some(tool.n_cpu),
        );
        let mut md = PileupMetadata::new();
        let mut start0: u32 = 0;
        for chunk in worker.pileup.0.chunk_by(|a, b| a == b){
            md.n_chunks += 1;
            let chunk_len = chunk.len();
            let end1 = start0 + chunk_len as u32;
            let chunk_coverage = chunk[0].coverage();
            if chunk_coverage > 0 {
                let excluded  =  tool.exclusions.span_overlaps_region(&worker.chrom, start0 + 1, end1);
                let on_target = !tool.targets.has_data || 
                                       tool.targets.span_overlaps_region(&worker.chrom, start0 + 1, end1);
                if !excluded && on_target {
                    let record = PileupRecord {
                        chrom_index: worker.chrom_index,
                        start0,
                        end1,
                        counts: &chunk[0],
                    };
                    csv.serialize(&record);
                    md.n_reported_chunks += 1;
                    md.n_reported_bases  += chunk_len;
                    md.reported_coverage += chunk_coverage * chunk_len;
                    if chunk_coverage > chunk[0].n_ref as usize {
                        md.n_variant_bases += chunk_len ;
                    }
                }
            }
            start0 = end1;
        }
        md
    }
}
