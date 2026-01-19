//! Tools and methods for grouping, i.e., fuzzy matching, junctions
//! into single SV calls.

// chr1
// +01000000001                            +01000000999
// --------------------------------------------------->
// <---------------------------------------------------
// -01000000001                            -01000000999
// chr2
// +02000000001                            +02000000999 >>> sorts last
// --------------------------------------------------->
// <---------------------------------------------------
// -02000000001                            -02000000999 <<< sorts first

// dependencies
use std::error::Error;
use rustc_hash::FxHashMap;
use rayon::prelude::*;
use genomex::{genome::Chroms, sam::SamRecord};
use super::{OrderedJunction, JunctionInstances, FinalJunction};

// constants
const LEFTWARD:  u8 = 0;
const RIGHTWARD: u8 = 1;

/// Write a vector of OrderedJunctions and associated JunctionInstances to a file, one row per
/// unique OrderedJunctions with concatenated fields and a count.
pub fn fuzzy_match_junctions(
    mut jxns_in: FxHashMap<OrderedJunction, JunctionInstances>,
    chroms:      &Chroms,
    group_breakpoint_distance: usize, // matching tolerances
    group_stem_distance:       u32,
) -> Result<Vec<FinalJunction>, Box<dyn Error>> {
    let mut jxns_out: Vec<FinalJunction> = Vec::new();

    // sort junctions based on first breakpoint node, ascending
    // bottom (-) strand sorts last chrom/pos to first, then top (+) strand sorts first chrom/pos to last
    let mut ordered_jxns: Vec<OrderedJunction> = jxns_in.keys().copied().collect();
    ordered_jxns.par_sort_unstable_by(|a, b| 
        a.node_aln5_end3.cmp(&b.node_aln5_end3)
    );

    // group junctions by first breakpoint node
    for node1_group in ordered_jxns.chunk_by_mut(|a, b| {
        (b.node_aln5_end3.abs() - a.node_aln5_end3.abs()).unsigned_abs() <= group_breakpoint_distance
    }){

        // if only one junction in group, write it as is
        if node1_group.len() == 1 {
            commit_junction(node1_group, &mut jxns_in, &mut jxns_out, chroms)?;
            continue;
        }

        // sort junctions based on second breakpoint node
        node1_group.sort_unstable_by(|a, b| {
            a.node_aln3_end5.cmp(&b.node_aln3_end5)
        });

        // group junctions by second breakpoint node
        for node2_group in node1_group.chunk_by_mut(|a, b| {
            (b.node_aln3_end5.abs() - a.node_aln3_end5.abs()).unsigned_abs() <= group_breakpoint_distance
        }){

            // if only one junction in group, write it as is
            if node2_group.len() == 1 {
                commit_junction(node2_group, &mut jxns_in, &mut jxns_out, chroms)?;
                continue;
            }

            // sort junctions based on var_chrom_len statistic
            node2_group.sort_unstable_by(|a, b| {
                var_chrom_len(a, chroms).cmp(&var_chrom_len(b, chroms))
            });

            // group junctions by var_chrom_len statistic
            for var_chrom_group in node2_group.chunk_by_mut(|a, b| {
                (var_chrom_len(b, chroms) - var_chrom_len(a, chroms)).unsigned_abs() <= group_stem_distance
            }){
                commit_junction(var_chrom_group, &mut jxns_in, &mut jxns_out, chroms)?;
            }
        }
    }
    Ok(jxns_out)
}

// commit one junction group from one or more fuzzy-matched OrderedJunctions
// even if only one OrderedJunction is in the group, it may have multiple JunctionInstances
// with different properties, e.g., different inserted bases, etc.
fn commit_junction(
    group_jxns: &[OrderedJunction],
    jxns_in:    &mut FxHashMap<OrderedJunction, JunctionInstances>,
    jxns_out:   &mut Vec<FinalJunction>,
    chroms:     &Chroms,
) -> Result<(), Box<dyn Error>> {

    // one junction in group, commit directly
    let final_jxn = if group_jxns.len() == 1 {
        FinalJunction::from_best_junction(
            &group_jxns[0], 
            jxns_in.get_mut(&group_jxns[0]).unwrap(),
            chroms
        )

    // pick the representative/best OrderedJunction as the one with the highest number of JunctionInstances
    // collect all instances starting from the best junction and commit
    } else {
        let best_jxn = group_jxns.iter().max_by_key(|jxn| {
            jxns_in[jxn].read_data.len()
        }).unwrap();
        let mut instances = jxns_in[best_jxn].clone();
        for jxn in group_jxns.iter().filter(|jxn| *jxn != best_jxn) {
            instances.extend(jxns_in.get_mut(jxn).unwrap());
        }
        FinalJunction::from_best_junction(
            &best_jxn, 
            &mut instances,
            chroms
        )
    };

    // store the resulting FinalJunction
    // TODO: ready to print instead?
    jxns_out.push(final_jxn);
    Ok(())
}

// estimate the number of genome bp still present as the chromosome stem length on each side of the junction
// 1     7                   9       1 (i.e., count rightward stems from the end of the chrom)
// -----------------------------------
//    1--> (7 bp stem)       <--1 (9 bp stem)
//    <--2 (7 bp stem)       2--> (9 bp stem)
// inappropriately clipped junction bases always do two things to ~the same degree:
//   decrease the summed chrom stem length
//   increase the size of the alignment offset, i.e., the number of apparently inserted bases
// thus, var_chrom_len = summed chrom stem length + alignment offset:
//   approximates the size of the resulting chromosome if this were the only junction present
//   stays ~constant as a clip-resistant metric, since most bases are accounted for in a read even if sequenced poorly
// var_chrom_len is a junction/SV-level property, unlike breakpoint distances which are assessed per breakpoint
// for efficiency, each breakpoint is assessed more promiscuously in series to establish candidate junction groups
// then the group is subjected to var_chrom_len matching more stringently since clip errors can now be accounted for
// it is not practical to account for clip errors during breakpoint distance matching
// because we cannot know which breakpoint (either, both, neither) may have been inappropriately clipped
// var_chrom_len, in considering both breakpoints at once, doesn't care which breakpoint was inappropriately clipped
// the resulting value is the same for:
//    no inappropriate clipping
//    read 1 inappropriately clipped
//    read 2 inappropriately clipped
//    both read 1 and read 2 inappropriately clipped
// demonstrating the resilience of this metric for accurate junction path matching
// one cannot use only var_chrom_len since 
//    it would be computationally inefficient,
//    wildly different junction positions could yield the same value
// but once breakpoints are known to be roughly co-localized, var_chrom_len provides robust junction matching
// the method is robust to inversions and all types of translocations
fn var_chrom_len(
    jxn:    &OrderedJunction,
    chroms: &Chroms,
) -> i32 {
    let (_, chrom_index_1, ref_pos1_1, is_reverse_1) = 
        SamRecord::unpack_signed_node(jxn.node_aln5_end3, chroms);
    let (_, chrom_index_2, ref_pos1_2, is_reverse_2) = 
        SamRecord::unpack_signed_node(jxn.node_aln3_end5, chroms);
    let chrom_stem_1 = chrom_stem_len(
        chrom_index_1, 
        ref_pos1_1,
        if is_reverse_1 { RIGHTWARD } else { LEFTWARD },
        chroms,
    );
    let chrom_stem_2 = chrom_stem_len(
        chrom_index_2, 
        ref_pos1_2,
        if is_reverse_2 { LEFTWARD } else { RIGHTWARD },
        chroms,
    );
    (chrom_stem_1 as i32 + chrom_stem_2 as i32 + jxn.offset as i32).abs() // should always be positive...
}
// get one chrom stem length
fn chrom_stem_len(
    chrom_index: u8,
    ref_pos1:    u32,
    side0:       u8,
    chroms:      &Chroms,
) -> u32 {
    if side0 == LEFTWARD {
        ref_pos1
    } else {
        chroms.index_sizes[&chrom_index].saturating_sub(ref_pos1) + 1
    }
}
