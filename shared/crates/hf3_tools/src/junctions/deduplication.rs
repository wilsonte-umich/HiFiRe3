//! Separate out the important logic of read deduplication for junction instances.

// dependencies
use rustc_hash::FxHashMap;
use super::{JunctionInstances, FinalJunction};

impl FinalJunction {

    /// Deduplicate ONT reads in a set of JunctionInstances by setting the 
    /// read_data.duplicate flag based on ONT channel and junction orientation.
    /// 
    /// ONT allows one read per channel per junction instance orientation
    /// regardless of outer endpoints. 
    /// 
    /// Opposite orientation reads in the same channel are nearly always duplex 
    /// sequencing of two strands from the same original DNA molecule and should 
    /// be deduplicated accordingly.
    /// 
    /// In contrast, two reads from the same channel in the same orientation
    /// sequenced the same DNA strand and therefore must have arisen from
    /// different original DNA molecules in a PCR-free library.
    pub (super) fn mark_duplicates_by_channel(
        instances: &mut JunctionInstances,
        cache:     &mut FxHashMap<u32, FxHashMap<u8, Vec<usize>>>,
    ){
        // cache = <read_data.channel, <instance.jxn_orientation, Vec<instance_index>>>
        cache.clear();
        for instance_i in 0..instances.read_data.len() {
            let channel = instances.read_data[instance_i].channel;
            let jxn_orientation = instances.jxn_orientations[instance_i];
            cache.entry(channel)
                .or_insert_with(FxHashMap::default)
                .entry(jxn_orientation)
                .or_insert_with(Vec::new)
                .push(instance_i);
        }
        // for each channel, choose the orientation with the most instances
        // mark all opposite-orientation detections as duplicates
        // it is possible that no unmarked instance will have offset == best_offset_abs
        for (_channel, ori_map) in cache.iter() {
            //   -------|---A0---|------> ch1    -------|---C0---|------> ch2
            //     <----|---A1---|-------                 <-C1---|-------     A1/C1 duplex reads are often incomplete, mark as dups regardless
            //      ----|---B0---|----->            ----|---D0---|----->      B0/D0 must be unique relative to A0/C0, do not mark as dups
            if ori_map.keys().len() == 1 { continue; } // nothing to deduplicate for this channel, only one strand present
            let best_ori = ori_map.keys().max_by_key(|ori| ori_map[ori].len()).unwrap();
            for instance_i in ori_map[&(1 - *best_ori)].iter() { // i.e., all instances in the opposite orientation
                instances.read_data[*instance_i].is_duplicate = true;
            }
        }
    }

    /// Deduplicate non-ONT, non-RE reads in a set of JunctionInstances by setting the 
    /// read_data.duplicate flag based on read outer endpoints, e.g., PCR-amplified  
    /// Nextera libraries allow only one recorded read per unique outer endpoint pair.
    ///  
    /// There is a general presumption that multiple reads with the same random outer 
    /// endpoints represent PCR duplicates of the same original DNA molecule, which 
    /// could arise on either strand. The number of cases where this presumption is 
    /// violated should be low, especially if the library complexity is high and the 
    /// sequencing depth is moderate, but will lead to undercounting when it occurs.
    /// 
    /// Reduced representation RE-mediated libraries are never deduplicated in this way 
    /// since outer endpoints are not random and multiple independent DNA molecules are 
    /// expected to share the same outer endpoints.
    /// 
    /// Could deploy fuzzy matching on outer endpoints to deduplicate more 
    /// aggressively, but short-read outer endpoints are mostly expected to be accurate.
    pub (super) fn mark_duplicates_by_span(
        instances: &mut JunctionInstances,
        cache:     &mut FxHashMap<(isize, isize), Vec<usize>>,
        best_offset_abs: usize,
    ){
        // cache = <(read_data.outer_node1, read_data.outer_node2), Vec<instance_index>>
        cache.clear();
        //    
        for instance_i in 0..instances.read_data.len() {
            let outer_endpoints = (
                instances.read_data[instance_i].node1, 
                instances.read_data[instance_i].node2
            );
            cache.entry(outer_endpoints)
                .or_insert_with(Vec::new)
                .push(instance_i);
        }
        for (_outer_endpoints, oe_is) in cache.iter(){
            //   -------|---A---|-------  two instances of the same junction in one read
            //   -------|---B---|-------  a second read with the same instances (the same or opposite orientation)
            if oe_is.len() == 1 { continue; } // nothing to deduplicate for this set of outer endpoints
            // keep the read from the first instance where offset.len() == best_offset_abs
            // if none, arbitrarily keep the first encountered qname
            let best_qname = oe_is.iter()
                .find(|&&i| instances.jxn_seqs[i].len() == best_offset_abs)
                .map(|&i|   instances.read_data[i].qname.clone())
                .unwrap_or_else(|| instances.read_data[oe_is[0]].qname.clone());
            // marks junction instances from all other reads as duplicates
            // multiple instances from the best QNAME all remain unmarked
            for instance_i in oe_is {
                if instances.read_data[*instance_i].qname != best_qname {
                    instances.read_data[*instance_i].is_duplicate = true;
                }
            }
        }
    }
}
