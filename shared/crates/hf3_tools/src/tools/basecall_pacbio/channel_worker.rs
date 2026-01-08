//! Handling strand merging, including decisions not to merge.
//! Generate outcomes for writing and tallying. Operates in parallel 
//! worker threads.

// modules
mod strand_merger;

// dependencies
use std::error::Error;
use crossbeam::{channel::{Receiver, Sender}};
use mdi::pub_key_constants;
use genomex::sequence::{rc_acgt_str, Aligner, AlignmentStatus};
use super::{StrandPair, KineticsInstance, MergeResult, MAX_READ_LEN};
use strand_merger::StrandMerger;

// constants
pub_key_constants!(
    FAILED_STRAND_ALIGNMENT
    ALIGNMENT_TOO_POOR
    ALIGNMENT_TOO_SHORT
    TWO_STRANDS_MERGED
);
const MAX_SHIFT: usize = 15; // allowed Smith-Waterman alignment shift from identity
const MAX_SCORE_DIFF_PER_KB: i32 = 25; // require a minimal match between strands ...
const MAX_SUMMED_CLIPS: usize = 25;    // ... that extends nearly to the read ends

// called by crossbeam scope to create parallel strand-merging worker threads
pub fn merge_strand_pairs(
    rx_strand_pair:  Receiver<StrandPair>,
    tx_kinetics:     Sender<KineticsInstance>,
    tx_merge_result: Sender<MergeResult>,
) -> Result<(), Box<dyn Error>> {

    // instantiate tools for this worker thread
    let mut aligner = Aligner::new_fast(
        MAX_READ_LEN + 1,
        MAX_READ_LEN + 1,
        MAX_SHIFT,
    );
    let mut strand_merger = StrandMerger::new();

    // process strand pairs as they arrive
    for strand_pair in rx_strand_pair.iter() {
        merge_strands(
            strand_pair, 
            &mut aligner, 
            &mut strand_merger,
            &tx_kinetics,
            &tx_merge_result,
        )?
    }
    Ok(())
}

// merge two usable strands into a single read, masking disagreements to Ns
fn merge_strands(
    strand_pair: StrandPair,
    aligner: &mut Aligner,
    strand_merger: &mut StrandMerger,
    tx_kinetics:     &Sender<KineticsInstance>,
    tx_merge_result: &Sender<MergeResult>,
) -> Result<(), Box<dyn Error>> {

    // reverse complement the previous strand sequence for alignment
    let qry = rc_acgt_str(&strand_pair.prev.seq);
    let tgt = &strand_pair.this.seq;

    // align the two strands
    let aln = aligner.align(&qry, tgt, None, true);

    // handle the unxpected outcome of failed alignment between strands
    if aln.status != AlignmentStatus::AlignmentFound {
        return abort_read(FAILED_STRAND_ALIGNMENT, strand_pair, &tx_merge_result);
    }

    // handle poor alignments between strands that might occur from failed strands or other poor ccs
    let max_len = tgt.len().max(qry.len()) as i32; // the expected perfect alignment score
    let min_allowed_score = max_len - (MAX_SCORE_DIFF_PER_KB * (max_len + 999) / 1000);
    if aln.score < min_allowed_score {
        return abort_read(ALIGNMENT_TOO_POOR, strand_pair, &tx_merge_result);
    }
    if aln.tgt_end0 - aln.tgt_start0 + 1 < max_len as usize - MAX_SUMMED_CLIPS {
        return abort_read(ALIGNMENT_TOO_SHORT, strand_pair, &tx_merge_result);
    }

    // merge the two strands into a single read
    let (seq, qual, dd_tag, sk_tag, dt_tag) = 
        strand_merger.merge_strands(&strand_pair, &aln, tx_kinetics);

    // transmit the merged read with the two-strand tag set
    tx_merge_result.send(MergeResult {
        outcome:  TWO_STRANDS_MERGED,
        qname:    strand_pair.qname,
        seq:      seq,
        qual:     qual,
        ff:       strand_pair.ff,
        ec:       strand_pair.this.ec + strand_pair.prev.ec,
        dd:       Some(dd_tag),
        sk:       Some(sk_tag),
        dt:       Some(dt_tag),
    })?;
    Ok(())
}

// transmit a single-read when merging is not possible
fn abort_read(
    outcome: &'static str,
    strand_pair: StrandPair,
    tx_merge_result: &Sender<MergeResult>
) -> Result<(), Box<dyn Error>> {
    tx_merge_result.send(MergeResult {
        outcome,
        qname:    strand_pair.qname,
        seq:      strand_pair.this.seq,
        qual:     strand_pair.this.qual,
        ff:       strand_pair.ff,
        ec:       strand_pair.this.ec,
        dd:       None,
        sk:       None,
        dt:       None,
    })?;
    Ok(())
}
