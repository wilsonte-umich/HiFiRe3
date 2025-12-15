//! This is a template for constructing a tool worker submodule using MDI framework components.
//! 
//! A tool worker is claimed as a dependency by a workflow tool and performs a specific
//! data processing task as part of the tool's overall workflow.
//! 
//! Replace this comment block with a description of the workers's purpose, actions
//! performed, expected inputs, and the generated outputs.

// dependencies
use mdi::pub_key_constants;
use mdi::workflow::Workflow;
use super::Tool;

// constants
pub_key_constants!{
    // from environment variables
    WORKER_OPTION 
    WORKER_FLAG
    // derived configuration values
    WORKER_DERIVED_OPTION
    // counter keys
    N_EVENTS
    N_BY_WORKER_KEY
    N_BY_WORKER_VALUE
}
const WORKER_CONSTANT: u8 = 1; // additional fixed values not exposed as options

//     EXPECTING_ENDPOINT_RE_SITES
//     REJECTING_JUNCTION_RE_SITES
//     // RE site matching
//     ENZYME_NAME
//     BLUNT_RE_TABLE
//     OVERHANG5_RE_TABLE
// // fillEnvVar(\our $SITE_CHROM_DATA_FILE,      'SITE_CHROM_DATA_FILE');     # in order of usage: first access the chrom's data
// // fillEnvVar(\our $CLOSEST_SITE_LOOKUP_WRK,   'CLOSEST_SITE_LOOKUP_WRK');  # then find the closest site on the chrom to query pos1
// // fillEnvVar(\our $SITE_DATA_LOOKUP_WRK,      'SITE_DATA_LOOKUP_WRK');     # then acquire the position and matching haplotypes
//     // tolerances and thresholds
//     CLIP_TOLERANCE      // RE site matching
//     ACCEPT_ENDPOINT_DISTANCE
//     REJECT_JUNCTION_DISTANCE

/// MySubModule does something useful.
pub struct SiteMatcher {
    my_value: bool,
}
impl SiteMatcher {
    /* ---------------------------------------------------------------------------
    initialize
    ---------------------------------------------------------------------------- */
    /// Initialize a new SiteMatcher.
    pub fn new(w: &mut Workflow) -> SiteMatcher {
        w.cfg.set_bool_env(&[WORKER_FLAG]);
        w.cfg.set_u8_env(&[WORKER_OPTION]);
        w.cfg.set_bool(WORKER_DERIVED_OPTION, *w.cfg.get_u8(WORKER_OPTION) == WORKER_CONSTANT);
        w.ctrs.add_counters(&[
            (N_EVENTS, "number of events processed"),
        ]);
        w.ctrs.add_keyed_counters(&[
            (N_BY_WORKER_KEY, "count of unique keys"),
        ]);
        w.ctrs.add_indexed_counters(&[
            (N_BY_WORKER_VALUE, "distribution of integer values"),
        ]);
        MySubModule{
            my_value: *w.cfg.get_bool(WORKER_FLAG),
        }
    }

    /* ---------------------------------------------------------------------------
    functions that do something useful
    ---------------------------------------------------------------------------- */
    /// Do something useful.
    pub fn do_something_useful(
        &self, 
        w: &mut Workflow, 
        tool: &mut Tool,
        record: &mut MyRecord,
    ) {

    }

}
