//!  smith_waterman() performs either Needleman-Wunsch or Smith-Waterman 
//!  alignment of a query (qry) to a target (tgt) sequence, depending on the arguments
//!  expects all upper case IUPAC codes
//!  returns just one best alignment as an array of query values
//!    M operations carry the query base in the array slot (could be a base mismatch)
//!    I operations carry the inserted base prepended to the NEXT target postion
//!    D operations carry "-" in place of the query base that was deleted relative to target
//!  in fast mode with max_shift > 0, query and target are assumed to be nearly identical and of similar length
//!    such as when comparing duplicate reads of the same sequence span
//!    thus, only a limited subset of all possible alignment registers are considered, set by max_shift
//!  in force_qry_end mode, alignment and score must always go to one (but not the other) end of query
//!    such as when placing a clip terminus into a candidate region on the other side of an SV junction
//!    the end of query to use is set by value of force_qry_end, either QRY_START or QRY_END
//!    this mode (and only this mode) will fail if there are multiple equally good alignments
//!  in local mode, best alignments can be local/partial, i.e., they do not need to go to the end of either molecule
//!  when not force_qry_end or local, alignments are inclusive to the end of at least one input sequence on each side
//!  with the defaults, max possible score is the length of the shorter input sequence
//!  because the expectation is that an Aligner will be reused for multiple alignments, the score matrix is pre-allocated
//!     based on max_expected_qry and max_expected_tgt provided at Aligner initialization
//!     be sure to use sufficiently large values to avoid panics due to memory access out of bounds
//!  for speed, callers can suppress creation of the base-by-base alignment map if only the score and end positions are needed 
//! 
//!  usage: 
//! ```
//! use genomex::sequence::smith_waterman::Aligner;
//! let mut aligner = Aligner::new(max_expected_qry, max_expected_tgt); // or with custom AlignerParameters via Aligner::with_param()
//! aligner.max_shift(5); // enable fast mode with max register shift of 5
//! aligner.suppress_alignment_map(); // suppress detailed alignment map generation
//! let alignment: Alignment = aligner.align(qry, tgt, Some(&ForceQryTerminus::QryStart), false); // force qry end alignment
//! if alignment.status != AlignmentStatus::AlignmentFound { ... }
//! println!("{}", alignment.score);
//! println!("{}", alignment.qry_on_tgt.join(""));
//! ```

// dependencies
use std::collections::HashMap;
use std::cmp::{min, max};
use super::iupac::*;

// constants for matrix traversal
const DIAG: u8 = 1;
const UP:   u8 = 2;
const LEFT: u8 = 3;

/// Enum to specify which end of the query to force alignment to, if any.
#[derive(PartialEq)]
pub enum ForceQryTerminus {
    QryStart,
    QryEnd,
}

/// Enum to specify alignment outcome status.
#[derive(PartialEq)]
pub enum AlignmentStatus {
    AlignmentFound,
    MultipleBestAlignments,
    NoValidAlignment,
}

/// Aligner scoring parameters.
pub struct AlignerParameters {
    match_score:            f32,   // first parameters have common defaults or can be set using with_param()
    mismatch_penalty:       f32,
    gap_open_penalty:       f32,
    gap_extension_penalty:  f32,
    max_shift:              usize, // next parameters have defaults but can be updated using max_shift() etc.
    suppress_alignment_map: bool,
    max_expected_qry:       usize, // last parameters are application-specific and must be set at Aligner initialization
    max_expected_tgt:       usize,
}

/// Alignment result.
pub struct Alignment {
    pub status:     AlignmentStatus,
    pub score:      i32,
    pub qry_start0: usize,
    pub qry_end0:   usize,
    pub tgt_start0: usize,
    pub tgt_end0:   usize,
    pub qry_on_tgt: Vec<String>, // will be empty if alignment map is suppressed
}

/// Smith-Waterman / Needleman-Wunsch aligner.
pub struct Aligner {
    pub param:       AlignerParameters,
    pub fast:        bool,
    pub pair_scores: HashMap<(u8, u8), f32>,
    pub matrix:      Vec<Vec<(f32, u8)>>, // pre-allocated score matrix
}
impl Aligner {

    /// Initialize an aligner with caller-provided AlignerParameters.
    pub fn with_param(param: AlignerParameters) -> Self {
        let fast = param.max_shift > 0;
        let pair_scores = Self::initalize_pair_scores(&param, Iupac::new());
        // score matrix of form matrix[target_index][query_index] = (score, pointer)
        // pre-allocate to max expected sizes of query and target sequences
        let matrix = vec![vec![(0.0, 0); param.max_expected_qry + 1]; param.max_expected_tgt + 1];
        Aligner {
            param,
            fast,
            pair_scores,
            matrix,
        }
    }

    /// Initialize an aligner with default AlignerParameters.
    pub fn new(max_expected_qry: usize, max_expected_tgt: usize) -> Self {
        Self::with_param(AlignerParameters {
            match_score:            1.0,
            mismatch_penalty:      -1.5,
            gap_open_penalty:      -2.501, // 0.001 ajustment gives slight preference to not opening a single-base terminal gap
            gap_extension_penalty: -1.0,
            max_shift:              0,     // no register shift limit by default
            suppress_alignment_map: false, // generate detailed alignment map by default, i.e., don't suppress it
            max_expected_qry,
            max_expected_tgt,
        })
    }

    /// Set the max_shift parameter to enable fast alignment mode when sequences are known to be in register.
    pub fn max_shift(&mut self, max_shift: usize) {
        self.param.max_shift = max_shift;
        self.fast = max_shift > 0;
    }

    /// Set whether to suppress generation of the detailed alignment map.
    pub fn suppress_alignment_map(&mut self) {
        self.param.suppress_alignment_map = true;
    }

    // create the lookup table for match/mismatch score for all possible IUPAC code combinations
    fn initalize_pair_scores(param: &AlignerParameters, iupac: Iupac) -> HashMap<(u8, u8), f32> {
        let mut scores: HashMap<(u8, u8), f32> = HashMap::new();
        iupac.base_matches.iter().for_each(|(&key, _)| {
            let chars: Vec<u8> = key.as_bytes().to_vec();
            scores.insert((chars[0], chars[1]), param.match_score);       // e.g. A:A, full match
        });
        iupac.ryswkm_matches.iter().for_each(|(&key, _)| {
            let chars: Vec<u8> = key.as_bytes().to_vec();
            scores.insert((chars[0], chars[1]), param.match_score / 2.0); // e.g. A:R, half match
        });
        iupac.n_matches.iter().for_each(|(&key, _)| {
            let chars: Vec<u8> = key.as_bytes().to_vec();
            scores.insert((chars[0], chars[1]), 0.0);                     // e.g. A:N, uninformative base position, neither promoted nor penalized
        });
        scores.insert((b'-', b'-'), 0.0);                                 // gap aligned to gap, no penalty
        let bases = vec!['A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'N', '-'];
        bases.iter().for_each(|&base1| {
            bases.iter().for_each(|&base2| {
                let key = (base1 as u8, base2 as u8);
                if !scores.contains_key(&key) {
                    scores.insert(key, param.mismatch_penalty);           // everything else defaults to a mismatch, with its penalty
                }
            });
        });
        scores
    }

    /// Run a Smith-Waterman / Needleman-Wunsch alignment.
    /// Expects all upper case IUPAC codes.
    pub fn align(
        &mut self,
        qry: &str, 
        tgt: &str, 
        force_qry_terminus: Option<&ForceQryTerminus>,
        local: bool
    ) -> Alignment {
        let make_map = !self.param.suppress_alignment_map;

        // check input sequences
        if qry.is_empty() { panic!("smith_waterman error: missing query/qry sequence"); }
        if tgt.is_empty() { panic!("smith_waterman error: missing target/tgt sequence"); }

        // shortcut identical sequences
        let n_q = qry.len();
        let n_t = tgt.len();
        if n_q == n_t && qry == tgt {
            return Alignment{
                status:     AlignmentStatus::AlignmentFound,
                score:      (self.param.match_score * n_q as f32) as i32,
                qry_start0: 0,
                qry_end0:   n_q - 1,
                tgt_start0: 0,
                tgt_end0:   n_t - 1,
                qry_on_tgt: if self.param.suppress_alignment_map { 
                    vec![]
                } else {
                    qry.chars().map(|c| c.to_string()).collect()
                }
            }
        }

        // prepare for alignment
        let mut forcing_to_qry_start = false;
        let is_forcing = if let Some(force_qry_terminus) = force_qry_terminus {
            forcing_to_qry_start = *force_qry_terminus == ForceQryTerminus::QryStart;
            true
        } else {
            false
        };
        let qry_rev: Vec<u8>;
        let tgt_rev: Vec<u8>;
        let (q, t) = if is_forcing && forcing_to_qry_start {
            // temporarily reverse sequence to allow same code for forceQryEnd QRY_START and QRY_END
            // algorithm is written for QRY_END (i.e., the right side of qry)
            qry_rev = qry.as_bytes().iter().rev().cloned().collect(); // must re-allocate
            tgt_rev = tgt.as_bytes().iter().rev().cloned().collect();
            (qry_rev.as_slice(), tgt_rev.as_slice())
        } else {
            (qry.as_bytes(), tgt.as_bytes())
        };
        for i_t in 0..=n_t { // zero out just the portions of the pre-allocated matrix to be used for this alignment
            self.matrix[i_t][0..=n_q].fill((0.0, 0));
        }
        let mut best_score: f32 = 0.0;
        let mut paths: Vec<(usize, usize)> = vec![];
        let mut best_path: (usize, usize) = (0, 0);

        // fill score matrix based on matches and gaps
        let (mut min_i_q, mut max_i_q) = (1, n_q);
        for i_t in 1..=n_t {
            let is_last_t = i_t == n_t;
            if self.fast { // limit the possible query to target base matches
                (min_i_q, max_i_q) = (
                    max(1,   i_t - self.param.max_shift), 
                    min(n_q, i_t + self.param.max_shift)
                )
            };
            for i_q in min_i_q..=max_i_q {
                let diag_score = self.matrix[i_t - 1][i_q - 1].0 + 
                    self.pair_scores[&(t[i_t - 1], q[i_q - 1])];
                let up_score = self.matrix[i_t - 1][i_q].0 +
                    if self.matrix[i_t - 1][i_q].1 == UP { // gap penalties
                        self.param.gap_extension_penalty
                    } else {
                        self.param.gap_open_penalty
                    };
                let left_score = self.matrix[i_t][i_q - 1].0 +
                    if self.matrix[i_t][i_q - 1].1 == LEFT {
                        self.param.gap_extension_penalty
                    } else {
                        self.param.gap_open_penalty
                    };
                let (score, pointer) = if diag_score >= up_score && diag_score >= left_score {
                    (diag_score, DIAG)
                } else if up_score >= diag_score && up_score >= left_score {
                    (up_score, UP)
                } else {
                    (left_score, LEFT)
                };
                if local || is_forcing {
                    self.matrix[i_t][i_q] = if score > 0.0 { (score, pointer) } else { (0.0, 0) }; // allow trimming of non forced query end
                    if local || i_q == n_q { // ensure that all reported alignments go to end of query, unless local
                        if score > best_score {
                            best_score = score;
                            paths = vec![(i_t, i_q)];
                        } else if score == best_score {
                            paths.push((i_t, i_q)); // array of equally good paths
                        }
                    }
                } else { // general untrimmed alignment when requiring end-to-end alignment, e.g. in consensus building
                    self.matrix[i_t][i_q] = (score, pointer); // best productive path to this residue pair
                    if (is_last_t || i_q == n_q) && score > best_score {
                        best_score = score;
                        best_path = (i_t, i_q); // just keep the first encountered path with the best score
                    }
                }


            }
        }

        // catch alignment failures due to multiple equally good alignments or no valid alignment
        if local || is_forcing { // demand just one best hit in force_qry_end mode
            if !local && paths.len() > 1 {
                return Alignment{
                    status:     AlignmentStatus::MultipleBestAlignments,
                    score:      0, // indicate failure due to multiple equally good alignments
                    qry_start0: 0,
                    qry_end0:   0,
                    tgt_start0: 0,
                    tgt_end0:   0,
                    qry_on_tgt: if make_map { vec!['!'.to_string(); n_q] } else { vec![] },
                }
            }
            best_path = paths[0]; // randomly take the first best path
        }
        if best_path.0 == 0 || best_path.1 == 0 {
            return Alignment{
                status:     AlignmentStatus::NoValidAlignment,
                score:      0, // indicate failure due to no valid alignment found
                qry_start0: 0,
                qry_end0:   0,
                tgt_start0: 0,
                tgt_end0:   0,
                qry_on_tgt: vec![],
            }
        }

        // trace backwards to deconvolute best matching path(s) and alignment map(s)
        let (mut i_t_max, mut i_q_max) = best_path;
        let mut qry_on_tgt: Vec<String> = vec![]; // will be built backwards, need to reverse below
        let (mut i_t, mut i_q) = (i_t_max, i_q_max);
        loop{
            let pointer = self.matrix[i_t][i_q].1;  
            if pointer == 0 { break; } // occurs either just after leftmost unaligned end or at beginning of a sequence
            if pointer == DIAG { // M operation relative to reference (could be a mismatched base)
                if make_map { qry_on_tgt.push((q[i_q - 1] as char).to_string()); }
                i_t -= 1;
                i_q -= 1;
            } else if pointer == LEFT { // I operation relative to reference, append to NEXT reference base
                if make_map {
                    if qry_on_tgt.is_empty() { // in case it is the first time through
                        qry_on_tgt.push((q[i_q - 1] as char).to_string());
                    } else {
                        qry_on_tgt.last_mut().map(|s| format!("{}{}", q[i_q - 1] as char, s));
                    };
                }
                i_q -= 1;
            } else {  // D operation relative to reference, pad with a dummy character 
                if make_map { qry_on_tgt.push("-".to_string()); }
                i_t -= 1;
            }
        }

        // revert back to original sequence orientation when QUERY_START
        if is_forcing && forcing_to_qry_start {
            i_q_max -= 1;
            i_t_max -= 1;
            i_q = n_q - i_q - 1;
            i_q_max = n_q - i_q_max - 1;
            i_t = n_t - i_t - 1;
            i_t_max = n_t - i_t_max - 1;
            (i_q, i_q_max, i_t, i_t_max) = (i_q_max, i_q, i_t_max, i_t);
            i_q_max += 1;
            i_t_max += 1;
        } else {
            qry_on_tgt.reverse();
        }

        // return best alignment
        Alignment{
            status:     AlignmentStatus::AlignmentFound,
            score:      best_score as i32,
            qry_start0: i_q, // 0-based positions of alignment endpoints on qry and tgt
            qry_end0:   i_q_max - 1,
            tgt_start0: i_t,
            tgt_end0:   i_t_max - 1,
            qry_on_tgt
        }
    }
}
