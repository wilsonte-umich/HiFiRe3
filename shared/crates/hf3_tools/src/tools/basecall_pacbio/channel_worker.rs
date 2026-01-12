//! Handling strand merging, including decisions not to merge.
//! Generate outcomes for writing and tallying. Operates in parallel 
//! worker threads.

// modules
mod strand_merger;

// dependencies
use std::error::Error;
use crossbeam::{channel::{Receiver, Sender}};
use minimap2::{Aligner as Minimap2, Built};
use faimm::IndexedFasta;
use mdi::pub_key_constants;
use super::{StrandPair, KineticsInstance, MergeResult, MAX_READ_LEN};
use strand_merger::StrandMerger;

// constants
pub_key_constants!(
    TWO_STRANDS_MERGED
);

// called by crossbeam scope to create parallel strand-merging worker threads
pub fn merge_strand_pairs(
    rx_strand_pair:  Receiver<StrandPair>,
    tx_kinetics:     Sender<KineticsInstance>,
    tx_merge_result: Sender<MergeResult>,
    minimap2:        &Minimap2<Built>,
    fa:              &IndexedFasta,
) -> Result<(), Box<dyn Error>> {

    // initialize tools for this worker thread
    let mut strand_merger = StrandMerger::new();

    // process strand pairs as they arrive
    for strand_pair in rx_strand_pair.iter() {

        // create prev_on_this map for discovering hetroduplex strands
        // abort and report just this strand if strands don't align sufficiently well
        if let Some(failure) = strand_merger.set_prev_on_this(&strand_pair) {
            tx_merge_result.send(MergeResult {
                outcome:  failure,
                qname:    strand_pair.qname,
                seq:      strand_pair.this.seq,
                qual:     strand_pair.this.qual,
                ff:       strand_pair.ff,
                ec:       strand_pair.this.ec,
                dd:       None,
                sk:       None,
                dt:       None,
            })?;
        } else {

            // create ref_on_this map for assigning is_ref to hetroduplex strands
            strand_merger.set_ref_on_this(&strand_pair, minimap2, fa);

            // merge the two strands into a single read
            let (seq, qual, dd_tag, sk_tag, dt_tag) = 
                strand_merger.merge_strands(&strand_pair, &tx_kinetics);

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
        }
    }
    Ok(())
}
