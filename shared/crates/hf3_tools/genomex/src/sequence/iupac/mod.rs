//! Tabulations of bases and base matches in IUPAC codes.
// A	Adenine
// C	Cytosine
// G	Guanine
// T (or U)	Thymine (or Uracil)
// R	A or G
// Y	C or T
// S	G or C
// W	A or T
// K	G or T
// M	A or C
// B	C or G or T
// D	A or G or T
// H	A or C or T
// V	A or C or G
// N	any base
// . or -	gap

// dependencies
use std::collections::HashMap;

// Data type with IUPAC base representation tables.
pub struct Iupac {
    pub base_matches:    HashMap<&'static str, char>,
    pub base_mismatches: HashMap<&'static str, char>,
    pub ryswkm_matches:  HashMap<&'static str, char>,
    pub n_matches:       HashMap<&'static str, char>,
    pub m_operations:    HashMap<&'static str, char>,
    pub acgtn_matches:   HashMap<&'static str, char>,
    pub no_indel_matches:HashMap<&'static str, char>,
    pub space_matches:   HashMap<&'static str, char>,
}

/* ------------------------------------------------------------------
declare tables of lookup values for all base combinations that have either:
    an output IUPAC code other than the default of N
    and output SW score other than the default of mismatchPenalty

    three-base degenerate codes are not used, only bases and 2-base degeneracy
    thus, a base set with three possible bases defaults to output N
    and code BDHV will never match anything

    e.g.  ACGTTCAARCC ----> ACRTTYAANCC
          ACATTYAAYCC    
------------------------------------------------------------------ */
impl Iupac {
    pub fn new() -> Self {
        // paired base values for incoming raw bases (A C G T)
        // pairs full matchScore in SW
        let base_matches = HashMap::from([
            ("AA", 'A'),
            ("CC", 'C'),
            ("GG", 'G'),
            ("TT", 'T'),
        ]);
        // pairs with full mismatchPenalty in SW
        let base_mismatches = HashMap::from([
            ("AG", 'R'),
            ("GA", 'R'),
            ("CT", 'Y'),
            ("TC", 'Y'),
            ("GC", 'S'),
            ("CG", 'S'),
            ("AT", 'W'),
            ("TA", 'W'),
            ("GT", 'K'),
            ("TG", 'K'),
            ("AC", 'M'),
            ("CA", 'M'),
        ]);
        // paired base values that include RYSWKM, i.e. partial degeneracy
        // pairs with half matchScore in SW
        let ryswkm_matches = HashMap::from([
            ("RR", 'R'),
            ("RA", 'R'),
            ("RG", 'R'),
            ("AR", 'R'),
            ("GR", 'R'),
            ("YY", 'Y'),
            ("YC", 'Y'),
            ("YT", 'Y'),
            ("CY", 'Y'),
            ("TY", 'Y'),
            ("SS", 'S'),
            ("SG", 'S'),
            ("SC", 'S'),
            ("GS", 'S'),
            ("CS", 'S'),
            ("WW", 'W'),
            ("WA", 'W'),
            ("WT", 'W'),
            ("AW", 'W'),
            ("TW", 'W'),
            ("KK", 'K'),
            ("KG", 'K'),
            ("KT", 'K'),
            ("GK", 'K'),
            ("TK", 'K'),
            ("MM", 'M'),
            ("MA", 'M'),
            ("MC", 'M'),
            ("AM", 'M'),
            ("CM", 'M'),
        ]);
        // paired base values for any combination that includes N, i.e. full degeneracy
        // score with neutral=0 in SW (neither penalize nor promote)
        // thus, a base previously declared uninformative is an M operation but has no alignment value
        let n_matches = HashMap::from([
            ("NN", 'N'),
            ("NA", 'N'),
            ("NC", 'N'),
            ("NG", 'N'),
            ("NT", 'N'),
            ("AN", 'N'),
            ("CN", 'N'),
            ("GN", 'N'),
            ("TN", 'N'),
            ("NR", 'N'),
            ("NY", 'N'),
            ("NS", 'N'),
            ("NW", 'N'),
            ("NK", 'N'),
            ("NM", 'N'),
            ("RN", 'N'),
            ("YN", 'N'),  
            ("SN", 'N'),
            ("WN", 'N'),
            ("KN", 'N'),
            ("MN", 'N'),
        ]); 
        // all combinations that do not default to output base N
        // used to generate consensus sequences
        let m_operations: HashMap<_, _> = base_matches.iter()
            .chain(base_mismatches.iter())
            .chain(ryswkm_matches.iter())
            .map(|(&k, &v)| (k, v)).collect();
        // combinations for doing ungapped comparisons via fn no_indel_match
        let acgtn_matches: HashMap<_, _> = base_matches.iter()
            .chain(n_matches.iter())
            .map(|(&k, &v)| (k, v)).collect();
        let no_indel_matches: HashMap<_, _> = base_matches.iter()
            .chain(base_mismatches.iter())
            .chain(n_matches.iter())
            .map(|(&k, &v)| (k, v)).collect();
        // combinations for assembling consensuses with missing data
        let space_matches = HashMap::from([
            (" A", 'A'),
            (" C", 'C'),
            (" G", 'G'),
            (" T", 'T'),
            (" N", 'N'),
            ("A ", 'A'),
            ("C ", 'C'),
            ("G ", 'G'),
            ("T ", 'T'),
            ("N ", 'N'),
            ("  ", ' '),
        ]);

        // assemble and return Iupac data structure
        Iupac{
            base_matches,
            base_mismatches,
            ryswkm_matches,
            n_matches,
            m_operations,
            acgtn_matches,
            no_indel_matches,
            space_matches,
        }
    }
}
