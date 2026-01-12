//! Codify the outermost nodes of reads to establish molecule identity.

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use genomex::sam::SamRecord;
use super::Tool;
use super::site_matching::SiteMatches;

// constants
pub_key_constants!{
    // from environment variables
    IS_COMPOSITE_GENOME
}

/// OuterNodes structure for determining read molecule identity.
pub struct OuterNodes {
    pub node5: isize,
    pub node3: isize,
}
impl OuterNodes {
    /// Determine the outer nodes for a set of alignments representing a read or read pair.
    pub fn get_outer_nodes(
        alns: &[SamRecord], 
        aln_sites: &mut [SiteMatches], 
        is_unmerged_pair: bool,
        paired_outer_node: Option<isize>,
        tool: &Tool,
        w: &mut Workflow,
    ) -> OuterNodes {
        let n_alns = alns.len();
        let (node5, node3): (isize, isize);

        // get outer nodes based on RE site (projections)
        if tool.expecting_endpoint_re_sites {
            // determine if the 5'-most alignment still needs to be projected
            if aln_sites[0].proj3.index1 == 0 {
                aln_sites[0].proj3 = tool.site_matcher.get_projection(&alns[0], aln_sites[0].site3.clone());
            }

            // determine if the 3'-most alignment still needs to be projected
            // here projection is more limited and simply bumps to the next RE site
            if aln_sites[n_alns - 1].proj3.index1 == 0 {
                aln_sites[n_alns - 1].proj3 = tool.site_matcher.get_projection(&alns[n_alns - 1], aln_sites[n_alns - 1].site3.clone());
            }

            // get outer nodes based on RE site positions
            // using site positions enables fuzzy endpoint matching when comparing reads
            node5 = alns[0].pack_signed_node_at_pos(aln_sites[0].site5.pos1, false, &tool.chroms); // reads always define their own 5' outer node
            node3 = if is_unmerged_pair {
                let (chrom, pos1, is_reverse) = SamRecord::unpack_signed_node(paired_outer_node.unwrap(), &tool.chroms);
                let closest_site = tool.site_matcher.find_closest_site(&chrom, pos1, w);
                SamRecord::pack_signed_node(&chrom, closest_site.pos1, !is_reverse, &tool.chroms) // invert the strand orientation of paired 5' ends to match the behavior of single or merged reads
            } else { // using projected 3' ends allows best deduplication of partially sequenced RE fragments
                alns[n_alns - 1].pack_signed_node_at_pos(aln_sites[n_alns - 1].proj3.pos1, false, &tool.chroms)
            };

        // get outer nodes based on alignment positions
        } else {
            node5 = alns[0].pack_signed_node_aln(5, &tool.chroms); // reads always define their own 5' outer node
            node3 = if is_unmerged_pair {
                let node = paired_outer_node.unwrap();
                -node // invert the strand orientation of paired 5' ends to match the behavior of single or merged reads
            } else { // this unprojectable non-RE node position may not be the end of the insert if a single read is not end-to-end
                alns[n_alns - 1].pack_signed_node_aln(3, &tool.chroms)
            };
        }
        OuterNodes { node5, node3 }
    }
    
    /// Encode outer nodes as a SAM tag string.
    pub fn to_sam_tag(&self) -> String {
        format!("{},{}", self.node5, self.node3)
    }
}
