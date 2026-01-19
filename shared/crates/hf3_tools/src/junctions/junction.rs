//! Structure and methods for aggregating and describing detected junctions.

// dependencies
// use std::io::Write;
use std::mem::take;
use rustc_hash::FxHashMap;
use rust_htslib::bam::record::Record as BamRecord;
use genomex::genome::Chroms;
use genomex::bam::tags;
use genomex::sam::{junction::Junction, SamRecord};
use crate::formats::hf3_tags::*;
use crate::inserts::{ReadLevelMetadata, UniqueInsertSpan};

/// An OrderedJunction describes a single detected junction in the canonical orientation.
/// 
/// Values are the same for all instances of a given junction and used to key (deduplicated) 
/// aggregation, where junction "sameness" is based on the breakpoint nodes and the number 
/// of bases in the junction microhomology or insertion. Thus, it considers the numbers
/// of bases in the junction path, but allows for base differences resulting from basecalling 
/// errors (but not indels).
#[derive(PartialEq, Eq, Hash, Clone, Copy)]
pub struct OrderedJunction {
    /* ------------------------------------------- */
    // two ordered breakpoint nodes and an offset define a junction
    // each of these values is needed for fuzzy matching and grouping of junctions
    pub node_aln5_end3: isize, // +|- chrom_index << 29 + junction node pos1
    pub node_aln3_end5: isize, // node_aln3_end5 = 5' end of the 3'-most alignment in canonical orientation
    pub offset:         i16, // + offset = insertion, - offset = microhomology overlap
    /* ------------------------------------------- */
    // junction types and sizes derive deterministically from nodes so are the same for all JunctionInstances
    pub jxn_type: u8,
    strands:      u8,
    sv_size:      u32, // modal size selected late when OrderedJunctions are aggregated by fuzzy matching
}

/// JunctionInstance is one observation of a given OrderedJunction from a read alignment.
pub struct JunctionInstance {
    /* ------------------------------------------- */
    // junction properties that collapse to a single value after junction aggregation
    /* ------------------------------------------- */
    jxn_seq: String, // bases in the canonical orientation; best value picked during but not used for junction grouping
    /* ------------------------------------------- */
    // values maintained as lists of values per junction instance
    /* ------------------------------------------- */
    // values that are the same for all read alignments (even those not flanking the junction)
    pub read_data: ReadLevelMetadata, // qname, insert_size, insert_span (node1, node2, channel), n_jxns
    /* ------------------------------------------- */
    // values that combine properties of the two flanking alignments
    jxn_orientation:  u8, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flag: u8, // JXN_FAILURE_FLAG tag value
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
#[derive(Clone)]
pub struct JunctionInstances {
    /* ------------------------------------------- */
    jxn_seqs: Vec<String>, // bases in the canonical orientation
    /* ------------------------------------------- */
    pub read_data: Vec<ReadLevelMetadata>, // qname, insert_size, insert_span (node1, node2, channel), n_jxns
    /* ------------------------------------------- */
    jxn_orientations:  Vec<u8>, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flags: Vec<u8>, // JXN_FAILURE_FLAG tag value
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
            jxn_seqs:          Vec::new(),
            /* ------------------------------------------- */
            read_data:         Vec::new(),
            /* ------------------------------------------- */
            jxn_orientations:  Vec::new(),
            jxn_failure_flags: Vec::new(),
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
        self.jxn_seqs         .push(instance.jxn_seq);
        /* ------------------------------------------- */
        self.read_data        .push(instance.read_data);
        /* ------------------------------------------- */
        self.jxn_orientations .push(instance.jxn_orientation);
        self.jxn_failure_flags.push(instance.jxn_failure_flag);
        self.aln_failure_flags.push(instance.aln_failure_flag);
        /* ------------------------------------------- */
        self.min_stem_lengths .push(instance.min_stem_length);
        self.min_mapqs        .push(instance.min_mapq);
        self.max_divergences  .push(instance.max_divergence);
    } 
    /// Add multiple JunctionInstances to the collection.
    pub fn extend(&mut self, other: &mut JunctionInstances) {
        /* ------------------------------------------- */
        self.jxn_seqs         .extend(take(&mut other.jxn_seqs));
        /* ------------------------------------------- */
        self.read_data        .extend(take(&mut other.read_data));
        /* ------------------------------------------- */
        self.jxn_orientations .extend(take(&mut other.jxn_orientations));
        self.jxn_failure_flags.extend(take(&mut other.jxn_failure_flags));
        self.aln_failure_flags.extend(take(&mut other.aln_failure_flags));
        /* ------------------------------------------- */
        self.min_stem_lengths .extend(take(&mut other.min_stem_lengths));
        self.min_mapqs        .extend(take(&mut other.min_mapqs));
        self.max_divergences  .extend(take(&mut other.max_divergences));
    }
}

/// Extract the OrderedJunction and JunctionInstance from the 5' alignment of a junction pair.
pub fn get_junction(
    aln5:      &BamRecord,
    aln3:      &BamRecord,
    read_data: &ReadLevelMetadata,
) -> (OrderedJunction, JunctionInstance) { // the 5' alignment carries the encoded junction metadata
    let jxn_tag = tags::get_tag_str(aln5, JUNCTION); // must exist if this function is called
    let jxn = Junction::deserialize(&jxn_tag).orient_junction();
    (
        OrderedJunction {
            /* ------------------------------------------- */
            node_aln5_end3: jxn.node_aln5_end3,
            node_aln3_end5: jxn.node_aln3_end5,
            offset:         jxn.offset,
            /* ------------------------------------------- */
            jxn_type:       jxn.jxn_type as u8,
            strands:        jxn.strands  as u8,
            sv_size:        jxn.sv_size  as u32,
        },
        JunctionInstance {
            /* ------------------------------------------- */
            jxn_seq:   jxn.jxn_seq,
            /* ------------------------------------------- */
            read_data: read_data.clone(),
            /* ------------------------------------------- */
            jxn_orientation:  jxn.was_flipped as u8,
            jxn_failure_flag: tags::get_tag_u8_default(aln5, JXN_FAILURE_FLAG, 0),
            aln_failure_flag: (
                tags::get_tag_u8_default(aln5, ALN_FAILURE_FLAG, 0) | 
                tags::get_tag_u8_default(aln3, ALN_FAILURE_FLAG, 0)
            ),
            /* ------------------------------------------- */
            min_stem_length:  tags::get_tag_u32_default(aln5, STEM_LENGTH5, 1).min(
                              tags::get_tag_u32_default(aln3, STEM_LENGTH3, 1)),
            min_mapq:         aln5.mapq().min(aln3.mapq()),
            max_divergence:   tags::get_tag_f32_default(aln5, DIVERGENCE, 0.0).max(
                              tags::get_tag_f32_default(aln3, DIVERGENCE, 0.0)),
        }
    )
}

/// FinalJunction describes on fully grouped and resolved junction.
pub struct FinalJunction {
    /* ------------------------------------------- */
    // two ordered breakpoint nodes and an offset define a junction
    chrom_index1_1:  u8,    // node1
    ref_pos1_1:      u32,
    strand_index0_1: u8,
    chrom_index1_2:  u8,    // node2
    ref_pos1_2:      u32,
    strand_index0_2: u8,
    offset:          i16, // junction offset
    jxn_seq:         String,
    /* ------------------------------------------- */
    // junction properties that derive deterministically from nodes
    jxn_type:        u8,
    strands:         u8,
    sv_size:         u32,
    /* ------------------------------------------- */
    n_observed:      u16, // number of independent junction instances observed
    /* ------------------------------------------- */
    // values maintained as lists of values per junction instance
    /* ------------------------------------------- */
    read_data: Vec<ReadLevelMetadata>, // qname, insert_size, insert_span (node1, node2, channel), n_jxns
    /* ------------------------------------------- */
    jxn_orientations:  Vec<u8>, // 0|1 bool if this junction instance was flipped relative the original read orientation
    jxn_failure_flags: Vec<u8>, // JXN_FAILURE_FLAG tag value
    aln_failure_flags: Vec<u8>, // bitwise OR of the aln-level flags on each side of the junction
    /* ------------------------------------------- */
    min_stem_lengths: Vec<u32>, // BEST stem length of a flanking alignment (only one has to pass)
    min_mapqs:        Vec<u8>,  // WORST mapping quality of a flanking alignment (a junction is only as good as its worst flank)
    max_divergences:  Vec<f32>, // WORST divergence of a flanking alignment
    /* ------------------------------------------- */
    // additional junction properties that derive deterministically from nodes, added later
    // UJXN_IS_INTERGENOME
    // UJXN_TARGET_1
    // UJXN_TARGET_DIST_1
    // UJXN_TARGET_2
    // UJXN_TARGET_DIST_2
    // UJXN_GENES_1
    // UJXN_GENE_DIST_1
    // UJXN_GENES_2
    // UJXN_GENE_DIST_2
    // UJXN_IS_EXCLUDED_1
    // UJXN_IS_EXCLUDED_2
    // UJXN_BKPT_COVERAGE_1
    // UJXN_BKPT_COVERAGE_2
    /* ------------------------------------------- */
    // additional read-level properties that dervie deterministacally from globacl config, added later
    // CMP_N_SAMPLES
    // CMP_SAMPLES    
}

impl FinalJunction {

    /// Convert a single best OrderedJunction and its JunctionInstances into a FinalJunction.
    /// Not all `instances` necessarily match exactly to `best_jxn` if they were fuzzy grouped.
    pub fn from_best_junction(
        best_jxn:  &OrderedJunction,
        instances: &mut JunctionInstances,
        chroms:    &Chroms,
    ) -> Self {

        // determine the best offset length for jxn_seq extraction
        let best_offset = best_jxn.offset.abs() as usize;

        // unpack the breakpoint nodes
        let (_, chrom_index1_1, ref_pos1_1, is_reverse_1) = 
            SamRecord::unpack_signed_node(best_jxn.node_aln5_end3, chroms);
        let (_, chrom_index1_2, ref_pos1_2, is_reverse_2) = 
            SamRecord::unpack_signed_node(best_jxn.node_aln3_end5, chroms);

        // deduplicate as needed
        // TODO: check the required deduplication logic, especially regarding ONT channels
        // ONT channels should allow only one instance of a junction per channer regardless of outer endpoints
        // Nextera libraries should allow only one instance based on outer endpoints regardless of channel
        if true { /////////////// TODO: need this condition check
            let mut by_span: FxHashMap<UniqueInsertSpan, usize> = FxHashMap::default();
            let mut to_remove: Vec<usize> = Vec::new();
            for i in 0..instances.read_data.len() {
                let span = &instances.read_data[i].insert_span;
                if let Some(j) = by_span.get(span).copied() {
                    if instances.jxn_seqs[j].len() != best_offset {
                        by_span.insert(*span, i);
                        to_remove.push(j);
                    } else {
                        to_remove.push(i);
                    }
                } else {
                    by_span.insert(*span, i);
                }
            }
            for i in to_remove.iter().rev() {
                instances.jxn_seqs         .swap_remove(*i);
                instances.read_data        .swap_remove(*i);
                instances.jxn_orientations .swap_remove(*i);
                instances.jxn_failure_flags.swap_remove(*i);
                instances.aln_failure_flags.swap_remove(*i);
                instances.min_stem_lengths .swap_remove(*i);
                instances.min_mapqs        .swap_remove(*i);
                instances.max_divergences  .swap_remove(*i);
            }
        }

        // get the most frequently observed jxn_seq value; must match the best_jxn offset
        let n_observed = instances.read_data.len() as u16;
        let best_jxn_seq = if n_observed == 1 {
            take(&mut instances.jxn_seqs[0]) // a singleton junction
        } else {
            instances.jxn_seqs.iter()
                .filter(|seq| seq.len() == best_offset)
                .fold(FxHashMap::default(), |mut acc, seq| {
                    *acc.entry(seq).or_insert(0) += 1;
                    acc
                })
                .into_iter()
                .max_by_key(|&(_, count)| count)
                .map(|(val, _)| val.as_str())
                .unwrap()
                .to_string()
        };

        // construct and return the FinalJunction
        FinalJunction {
            /* ------------------------------------------- */
            chrom_index1_1,
            ref_pos1_1,
            strand_index0_1: is_reverse_1 as u8,
            chrom_index1_2,
            ref_pos1_2,
            strand_index0_2: is_reverse_2 as u8,
            offset:          best_jxn.offset,
            jxn_seq:         best_jxn_seq,
            /* ------------------------------------------- */
            jxn_type:        best_jxn.jxn_type,
            strands:         best_jxn.strands,
            sv_size:         best_jxn.sv_size,
            /* ------------------------------------------- */
            n_observed,
            /* ------------------------------------------- */
            read_data:         take(&mut instances.read_data),
            /* ------------------------------------------- */
            jxn_orientations:  take(&mut instances.jxn_orientations),
            jxn_failure_flags: take(&mut instances.jxn_failure_flags),
            aln_failure_flags: take(&mut instances.aln_failure_flags),
            /* ------------------------------------------- */
            min_stem_lengths:  take(&mut instances.min_stem_lengths),
            min_mapqs:         take(&mut instances.min_mapqs),
            max_divergences:   take(&mut instances.max_divergences),
        }
    }

    // /// Write a vector of OrderedJunctions and associated JunctionInstances to a file, one row per
    // /// unique OrderedJunctions with concatenated fields and a count.
    // pub fn write_sorted_with_count(
    //     jxns: FxHashMap<OrderedJunction, JunctionInstances>,
    //     file_path: &str
    // ) -> Result<(), Box<dyn Error>> {
    //     let mut bgz_writer = BgzWriter::from_path(file_path)?;

    //     // sort junctions based on first breakpoint
    //     let mut ordered_jxns = jxns.keys().collect::<Vec<&OrderedJunction>>();
    //     ordered_jxns.par_sort_unstable_by(|a, b| 
    //         a.node_aln5_end3.cmp(&b.node_aln5_end3)
    //     );
    //     // let mut level1_indices: Vec<usize> = (0..ordered_jxns.len()).collect();

    //     // fuzzy group junctions by first breakpoint
    //     for node1_group in ordered_jxns.chunk_by_mut(|a, b| {
    //         b.node_aln5_end3 - a.node_aln5_end3 <= 10 // TODO: implement fuzzy matching threshold
    //     }){

    //         // if only one junction in group, write it as is

    //         // sort junctions based on second breakpoint
    //         node1_group.sort_unstable_by(|a, b| {
    //             a.node_aln3_end5.cmp(&b.node_aln3_end5)
    //         });

    //         // fuzzy group junctions by second breakpoint
    //         for node2_group in node1_group.chunk_by_mut(|a, b| {
    //             b.node_aln3_end5 - a.node_aln3_end5 <= 10 // TODO: implement fuzzy matching threshold
    //         }){

    //             // if only one junction in group, write it as is

    //             // enforce third grouping level based on stem whatever that was

    //         }

    //     }


    //     Ok(())
    // }   

    // // write a set of one or more junctions instances
    // fn write_junction(
    //     jxns: FxHashMap<OrderedJunction, JunctionInstances>,
    //     group_jxns: &[OrderedJunction], // may have one or more junctions, each with one or more instances
    //     bgz_writer:  &mut BgzWriter,
    // ) -> Result<(), Box<dyn Error>> {

    //     // first determine which OrderedJunction is the representative junction

    //     // then aggregate the JunctionInstances across all OrderedJunctions in the group

    //     // for ordered_jxn in ordered_jxns {

    //         // writeln!(
    //         //     bgz_writer, 
    //         //     "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
    //         //     seg.chrom_index1,
    //         //     seg.ref_pos5_observed,
    //         //     seg.ref_pos3_observed,
    //         //     seg.ref_pos3_projected,
    //         //     seg.strand_index0,
    //         //     seg.jxn_types,
    //         //     seg.n_jxns,
    //         //     seg.aln_i,
    //         //     chunk.len()
    //         // )?;        
    //     // }    
    //     Ok(())
    // }
}
