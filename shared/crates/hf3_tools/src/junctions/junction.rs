//! Structure and methods for aggregating and describing detected junctions.

// dependencies
use std::error::Error;
use std::io::Write;
use rust_htslib::bgzf::Writer as BgzWriter;
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bam::record::Record as BamRecord;
use rayon::prelude::*;
use genomex::bam::tags;
use genomex::sam::junction::{Junction, OrderedNodes};
use crate::formats::hf3_tags::*;

// constants

/// An OrderedJunction describes a single detected junction in the canonical orientation.
/// Values are the same for all instances of a given junction and used to key (deduplicated) aggregation.
#[derive(PartialEq, Eq, Hash)]
pub struct OrderedJunction {

    // two ordered breakpoint nodes define a junction
    node_aln5_end3: isize, // +|- chrom_index << 9 + junction node pos1
    node_aln3_end5: isize, // node_aln3_end5 = 5' end of the 3'-most alignment in canonical orientation

    // junction types and sizes derive deterministically from nodes so are the same for all JunctionInstances
    pub jxn_type: u8,
    strands:  u8,
    sv_size:  usize, // modal size will be selected when OrderedJunctions are aggregated by fuzzy matching

    // junctions detected in different channels are considered distinct for deduplication
    channel: u32,
}

/// JunctionInstance is one observation of a given OrderedJunction in a read alignment.
pub struct JunctionInstance {
    /* ------------------------------------------- */
    // values that collapse to a single value after junction aggregation
    offset:  isize,  // + offset = insertion, - offset = microhomology overlap
    jxn_seq: String, // in the canonical orientation
    /* ------------------------------------------- */
    // values maintained as list of values per junction instance
    jxn_orientation:  u8, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flag: u8, // JXN_FAILURE_FLAG tag value
    /* ------------------------------------------- */
    // values that are the same for all read alignments (even those not flanking the junction)
    pub qname:       String,
    pub insert_size: i32, // INSERT_SIZE tag value
    outer_nodes: OrderedNodes, // reoriented signed outer nodes used for deduplication; reflect site_pos1 when applicable
    n_jxns:      u8,  // number of junctions in the read path containing this alignment segment
    /* ------------------------------------------- */
    // values that aggregate the two flanking alignments into a single metadata metric
    /* ------------------------------------------- */
    // values that combine the values from the two flanking alignments
    aln_failure_flag: u8, // bitwise OR of the aln-level flags on each side of the junction
    /* ------------------------------------------- */
    // values that select one value from the two flanking alignments
    min_stem_length: u32, // BEST stem length of a flanking alignment (only one has to pass)
    min_mapq:        u8,  // WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergence:  f32, // WORST divergence of a flanking alignment
}

/// JunctionInstances collects lists of individual observations of a given OrderedJunction.
/// Values might differ between instances of the same junction, including offset and jxn_seq
/// which can depend on the sequenced insert bases, not just the breakpoint nodes.
pub struct JunctionInstances {
    /* ------------------------------------------- */
    offsets:  Vec<isize>,  // + offset = insertion, - offset = microhomology overlap
    jxn_seqs: Vec<String>, // in the canonical orientation
    /* ------------------------------------------- */
    jxn_orientations:  Vec<u8>, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flags: Vec<u8>, // JXN_FAILURE_FLAG tag value
    /* ------------------------------------------- */
    qnames:       Vec<String>,
    insert_sizes: Vec<i32>, // INSERT_SIZE tag value
    outer_nodes:  Vec<OrderedNodes>, // reoriented signed outer nodes used for deduplication; reflect site_pos1 when applicable
    n_jxns:       Vec<u8>,  // number of junctions in the read path containing this alignment segment
    /* ------------------------------------------- */
    aln_failure_flags: Vec<u8>, // bitwise OR of the aln-level flags on each side of the junction
    /* ------------------------------------------- */
    min_stem_lengths: Vec<u32>, // BEST stem length of a flanking alignment (only one has to pass)
    min_mapqs:        Vec<u8>,  // WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergences:  Vec<f32>, // WORST divergence of a flanking alignment
}
impl JunctionInstances {
    /// Create a new empty set of JunctionInstances. 
    pub fn new() -> Self {
        // do not pre-allocate vectors since number of instances is unknown and typically small
        JunctionInstances {
            /* ------------------------------------------- */
            offsets:           Vec::new(),
            jxn_seqs:          Vec::new(),
            /* ------------------------------------------- */
            jxn_orientations:  Vec::new(),
            jxn_failure_flags: Vec::new(),
            /* ------------------------------------------- */
            qnames:            Vec::new(),
            insert_sizes:      Vec::new(),
            outer_nodes:       Vec::new(),
            n_jxns:            Vec::new(),
            /* ------------------------------------------- */
            aln_failure_flags: Vec::new(),
            /* ------------------------------------------- */
            min_stem_lengths:  Vec::new(),
            min_mapqs:         Vec::new(),
            max_divergences:   Vec::new(),
        }
    }
    /// Add a JunctionInstance to the collection.
    pub fn add_instance(&mut self, instance: JunctionInstance) {
        /* ------------------------------------------- */
        self.offsets          .push(instance.offset);
        self.jxn_seqs         .push(instance.jxn_seq);
        /* ------------------------------------------- */
        self.jxn_orientations .push(instance.jxn_orientation);
        self.jxn_failure_flags.push(instance.jxn_failure_flag);
        /* ------------------------------------------- */
        self.qnames           .push(instance.qname);
        self.insert_sizes     .push(instance.insert_size);
        self.outer_nodes      .push(instance.outer_nodes);
        self.n_jxns           .push(instance.n_jxns);
        /* ------------------------------------------- */
        self.aln_failure_flags.push(instance.aln_failure_flag);
        /* ------------------------------------------- */
        self.min_stem_lengths .push(instance.min_stem_length);
        self.min_mapqs        .push(instance.min_mapq);
        self.max_divergences  .push(instance.max_divergence);
    } 
}

/// Extract the OrderedJunction and JunctionInstance from the 5' alignment of a junction pair.
pub fn get_junction(
    aln5:          &BamRecord,
    aln3:          &BamRecord,
    n_alns:        usize,
    ordered_outer_nodes: &OrderedNodes,
) -> (OrderedJunction, JunctionInstance) { // the 5' alignment carries the encoded junction metadata
    let channel    = tags::get_tag_u32_default(aln5, CHANNEL, 0);
    let jxn_tag = tags::get_tag_str(aln5, JUNCTION); // must exist if this function is called
    let jxn = Junction::deserialize(&jxn_tag).orient_junction();
    unsafe { (
        OrderedJunction {
            node_aln5_end3: jxn.node_aln5_end3,
            node_aln3_end5: jxn.node_aln3_end5,
            jxn_type:       jxn.jxn_type as u8,
            strands:        jxn.strands as u8,
            sv_size:        jxn.sv_size,
            channel,
        },
        JunctionInstance {
            offset:           jxn.offset,
            jxn_seq:          jxn.jxn_seq,
            /* ------------------------------------------- */
            jxn_orientation:  jxn.was_flipped as u8,
            jxn_failure_flag: tags::get_tag_u8_default(aln5, JXN_FAILURE_FLAG, 0),
            /* ------------------------------------------- */
            qname:            from_utf8_unchecked(aln5.qname()).to_string(),
            insert_size:      tags::get_tag_i32_default(aln5, INSERT_SIZE, aln5.seq_len() as i32),
            outer_nodes:      *ordered_outer_nodes,
            n_jxns:           n_alns as u8 - 1,
            /* ------------------------------------------- */
            aln_failure_flag: (
                tags::get_tag_u8_default(aln5, ALN_FAILURE_FLAG, 0) | 
                tags::get_tag_u8_default(aln3, ALN_FAILURE_FLAG, 0)
            ),
            /* ------------------------------------------- */
            min_stem_length:  tags::get_tag_u32(aln5, STEM_LENGTH5).min(
                              tags::get_tag_u32(aln3, STEM_LENGTH3)),
            min_mapq:         aln5.mapq().min(aln3.mapq()),
            max_divergence:   tags::get_tag_f32_default(aln5, DIVERGENCE, 0.0).max(
                              tags::get_tag_f32_default(aln3, DIVERGENCE, 0.0)),
        }
    ) }
}

/// Write a vector of OrderedJunctions and associated JunctionInstances to a file, one row per
/// unique OrderedJunctions with concatenated fields and a count.
pub fn write_sorted_with_count(
    jxns: FxHashMap<OrderedJunction, JunctionInstances>,
    file_path: &str
) -> Result<(), Box<dyn Error>> {
    let mut bgz_writer = BgzWriter::from_path(file_path)?;

    // perform initial sort of junctions
    let mut ordered_jxns = jxns.keys().collect::<Vec<&OrderedJunction>>();
    ordered_jxns.par_sort_unstable_by(|a, b| 
        a.node_aln5_end3.cmp(&b.node_aln5_end3)
    );
    let mut node1_indices: Vec<usize> = (0..ordered_jxns.len()).collect();

    // fuzzy group junctions by first breakpoint
    for node2_indices in node1_indices.chunk_by_mut(|a, b| a == b) { // TODO: implement fuzzy matching

        // fuzzy group junctions by second breakpoint
        node2_indices.sort_unstable_by(|a, b| {
            ordered_jxns[*a].node_aln3_end5.cmp(&ordered_jxns[*b].node_aln3_end5)
        });

    }

    for ordered_jxn in ordered_jxns {

        // writeln!(
        //     bgz_writer, 
        //     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
        //     seg.chrom_index1,
        //     seg.ref_pos5_observed,
        //     seg.ref_pos3_observed,
        //     seg.ref_pos3_projected,
        //     seg.strand_index0,
        //     seg.jxn_types,
        //     seg.n_jxns,
        //     seg.aln_i,
        //     chunk.len()
        // )?;        
    }
    Ok(())
}   

// write a set of one or more junctions instances
fn write_junction(
    jxns: FxHashMap<OrderedJunction, JunctionInstances>,
    ordered_jxn: &OrderedJunction,
    instances:   &JunctionInstances,
    bgz_writer:  &mut BgzWriter,
) -> Result<(), Box<dyn Error>> {
    Ok(())
}

