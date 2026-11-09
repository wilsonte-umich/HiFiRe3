//! Annotate SVs in a Sniffles2-compatible VCF with the size of the RE fragments 
//! that contain the junction for a requested series of REs, as well as the 
//! distance from the junction to the nearest RE site. Steps are to:
//!     - extract SV lines from the VCF
//!     - reject SVs that are too small or too close to chrom ends
//!     - collect all unique QNAMEs that supported any passing SV
//!     - collect SEQ for those reads from the BAM file that generated the VCF
//!     - for each SV, assemble its reference sequence, i.e., scaffold, from VCF data
//!     - use mini_consensus to assemble a read consensus for the SV allele, i.e., haplotype
//!     - locate the junction(s) in the consensus, noting that INV SVs have two junctions
//!     - in silico digest the consensus with each enzyme and report for one or two SV junctions:
//!         - any de novo base insertions in INFO/SVINSSEQS
//!         - fragment sizes in INFO/<RE_NAME>SIZES
//!         - nearest RE site distances in INFO/<RE_NAME>DISTS

// imports
use std::error::Error;
use std::str::from_utf8_unchecked;
use rustc_hash::FxHashMap;
use rust_htslib::bcf;
use rust_htslib::bcf::Read as VcfRead;
use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;
use rust_htslib::faidx;
use minimap2::{Aligner as Minimap2, Mapping};
use mdi::pub_key_constants;
use mdi::workflow::{Workflow, Config, Counters};
use mini_consensus::{Resolver, ResolverConfig};
use crate::benchmark::*;

// constants for environment variable, config, and counter keys, etc.
const TOOL: &str = "digest_sv_alleles";
pub_key_constants!(
    // from environment variables
    N_CPU
    GENOME_FASTA
    BENCHMARK_BAM
    SNIFFLES_VCF
    RESTRICTION_ENZYMES
    MIN_SV_LEN
    // counter keys
    N_INPUT_SVS
    N_OUTPUT_SVS
    N_UNIQUE_QNAMES
);
/// Tool for processing VCF SV calls.
struct DigestionTool {
    n_cpu: usize,
    min_sv_len: u32,
    svs:  Vec<bcf::Record>,
    scfs: Vec<Vec<Vec<u8>>>,
    seqs: FxHashMap<Vec<u8>, Vec<u8>>,
    res:  Vec<ReSite>,
}

// function called by hf3_tools main()
pub fn main() -> Result<(), Box<dyn Error>> {

    // get config from environment variables
    let mut cfg = Config::new();
    cfg.set_string_env(&[
        GENOME_FASTA, BENCHMARK_BAM, SNIFFLES_VCF, 
        RESTRICTION_ENZYMES
    ]);
    cfg.set_usize_env( &[N_CPU, MIN_SV_LEN]);

    // initialize counters
    let ctrs = Counters::new(TOOL, &[
        (N_INPUT_SVS,     "number of input SV lines in VCF"),
        (N_OUTPUT_SVS,    "number of output SV lines after size filtering"),
        (N_UNIQUE_QNAMES, "number of unique reads used to assemble SV consensuses"),
    ]);

    // initialize the tool
    let mut w = Workflow::new(TOOL, cfg, ctrs);
    w.log.initializing();

    // create the digestion tool
    let mut tool = DigestionTool {
        n_cpu:      *w.cfg.get_usize(N_CPU),
        min_sv_len: *w.cfg.get_usize(MIN_SV_LEN) as u32,
        svs:  Vec::with_capacity(10_000),
        scfs: Vec::with_capacity(10_000),
        seqs: FxHashMap::default(),
        res: w.cfg.get_string(RESTRICTION_ENZYMES).split(',').map(|re: &str|{
            let re: Vec<_> = re.split(':').collect();
            ReSite::from_re_name(re[0], Some(re[1]))
        }).collect(),
    };
    
    // execute the workflow
    let mut header = scan_vcf(&mut w, &mut tool)?;
    collect_read_seqs(&mut w, &mut tool)?;
    let mut vcf = get_bcf_writer(&tool, &mut header)?;
    resolve_consensuses(&mut w, &mut tool, &mut vcf)?;

    // print counts
    w.ctrs.print_grouped(&[
        &[N_INPUT_SVS, N_OUTPUT_SVS, N_UNIQUE_QNAMES],
    ]);
    Ok(())
}

// -----------------------------------------------------------------------------
// initial VCF scan
// -----------------------------------------------------------------------------
/// Scan the input VCF file to collect SV calls and unique QNAMEs.
fn scan_vcf (w: &mut Workflow, tool: &mut DigestionTool) -> Result<bcf::Header, Box<dyn Error>> {
    w.log.print("scanning input VCF");

    // initialize all readers
    let sniffles_vcf = w.cfg.get_string(SNIFFLES_VCF).to_string();
    let mut vcf = bcf::Reader::from_path(sniffles_vcf)?;
    let header = vcf.header().to_owned();
    let genome_fasta = w.cfg.get_string(GENOME_FASTA).to_string();
    let faidx = faidx::Reader::from_path(genome_fasta)?;

    // process SV records one at a time
    // ignore any non-SV variants, although they are not expected from benchmark VCF files
    for result in vcf.records() {
        let sv = result?;
        let Some(sv_type) = get_vcf_sv_type(&sv)? else { continue };
        w.ctrs.increment(N_INPUT_SVS);

        // skip SVs that are too small
        if sv_type != b"BND" {
            let Some(sv_len) = get_vcf_sv_len(&sv)? else { continue };
            let min_sv_len: u32 = if sv_type == b"INV" {
                SCAFFOLD_FLANK_LEN as u32 // ensure that scaffolds don't cross both INV junctions
            } else {
                tool.min_sv_len
            };
            if sv_len < min_sv_len { continue }     
        }

        // collect information required to collect SEQs and build consensus for site analysis
        let Some(qnames) = get_vcf_qnames(&sv)? else { continue };
        if let Some(scfs) = build_scaffolds(&faidx, &sv, &sv_type)? {
            tool.svs.push(sv);   
            tool.scfs.push(scfs); 
            let downsampled_qnames = &qnames[0..qnames.len().min(MAX_CONSENSUS_SEQS)];
            for qname in downsampled_qnames {
                tool.seqs.insert(qname.clone(), Vec::new());
            }       
        };
    };

    // return the input header for modifcation by the VCF writer
    w.ctrs.add_to(N_OUTPUT_SVS, tool.svs.len());
    w.ctrs.add_to(N_UNIQUE_QNAMES, tool.seqs.len());
    Ok(bcf::Header::from_template(&header))
}

/// Build the reference scaffold sequence(s) for an SV's junction(s).
fn build_scaffolds(
    // header:  &bcf::header::HeaderView,
    faidx:   &faidx::Reader, 
    sv:      &bcf::Record,
    sv_type: &[u8],
) -> Result<Option<Vec<Vec<u8>>>, Box<dyn Error>> {
    let alt  = get_vcf_alt(&sv)?;

    // collect metadata on the left side of the SV
    let (chrom_1, pos1_1) = get_vcf_left_node(&sv)?;
    let chrom_len_1 = faidx.fetch_seq_len(&chrom_1) as usize;
    if pos1_1 <= SCAFFOLD_FLANK_LEN || pos1_1 + SCAFFOLD_FLANK_LEN >= chrom_len_1 { 
        return Ok(None) 
    }

    // collect metadata on the right side of the SV
    let (chrom_2, pos1_2) = get_vcf_right_node(&sv, &alt, &chrom_1)?;
    let chrom_len_2 = faidx.fetch_seq_len(&chrom_2) as usize;
    if pos1_2 <= SCAFFOLD_FLANK_LEN || pos1_2 + SCAFFOLD_FLANK_LEN >= chrom_len_2 { 
        return Ok(None) 
    }

    // parse the scaffold sequences for later use
    //  - throughout, Sniffles2 VCF does not provide information on bases inserted 
    //    at junctions, but these bases should resolve in the consensus
    //  - INFO/STRAND is NOT information about the junction path, it summarizes 
    //    the orientation(s) of supporting reads    
    let mut scfs: Vec<Vec<u8>> = Vec::new();
    match sv_type {
        b"DEL" => {
            // POS: 2707, REF: TCTgca...ccg, ALT: T, SVLEN=-259; END=2966; (2707 + 259 = 2966)
            // DEL is a single junction that is a fully resolved SV (sub)path
            let mut scf_left = get_ref_seq(faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
            let scf_right = get_ref_seq(faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;
            scf_left.extend(scf_right);
            scfs.push(scf_left);
        },
        b"INS" => {
            // POS: 2566, REF: c, ALT: cCCTCACCCTCACCCTCACCCTCACCCTGACCCTGACCCTGACCCTCACCCTGACCCTCACCCTCACCCTCA, SVLEN=71; END=2566;
            // INS is a single junction that is a fully resolved SV (sub)path (but DEL+INS at same junction is ambiguous)
            let mut scf_left = get_ref_seq(faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
            let scf_ins = alt[1..].to_ascii_uppercase();
            let scf_right = get_ref_seq(faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;
            scf_left.extend(scf_ins);
            scf_left.extend(scf_right);
            scfs.push(scf_left);
        },
        b"DUP" => {
            // POS: 3009031, REF: G, ALT: <DUP>, SVLEN=1394; END=3010425; (3009031 + 1394 = 3010425)
            // DUP is a single junction that is a fully resolved SV (sub)path
            let mut scf_left = get_ref_seq(faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, false)?;
            let scf_right = get_ref_seq(faidx, &chrom_1, pos1_1 + 1, pos1_1 + SCAFFOLD_FLANK_LEN, false)?;
            scf_left.extend(scf_right);
            scfs.push(scf_left);
        },
        b"INV" => {
            // POS: 40653961, REF: c, ALT: <INV>, SVLEN=434; END=40654395; (40653961 + 434 = 40654395)
            // INV is TWO junctions from a fully resolved inversion segment 

            // left junction
            let mut scf_left = get_ref_seq(faidx, &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
            let scf_right = get_ref_seq(faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, true)?;
            scf_left.extend(scf_right);
            scfs.push(scf_left);

            // right junction
            let mut scf_left = get_ref_seq(faidx, &chrom_1, pos1_1 + 1, pos1_1 + SCAFFOLD_FLANK_LEN, true)?;
            let scf_right = get_ref_seq(faidx, &chrom_2, pos1_2 + 1, pos1_2 + SCAFFOLD_FLANK_LEN, false)?;
            scf_left.extend(scf_right);
            scfs.push(scf_left);
        },
        b"BND" => {
            // CHR1: chr5, POS1: 16396441, REF: C, ALT: [chr1:145169227[C, CHR2=chr1; (no SVLEN or END)
            // orientation encoded in ALT
            // BND is a single junction from either a fully or partially resolved SV (sub)path

            // REF  ALT   Meaning
            // s    t[p[  piece extending to the right of p is joined after t
            // s    t]p]  reverse comp piece extending left of p is joined after t    *****
            // s    ]p]t  piece extending to the left of p is joined before t
            // s    [p[t  reverse comp piece extending right of p is joined before t  *****

            // I  do not know why mostly the ***** orientation are found in Sniffles2 output
            // found exactly one exception: chr10   t       ]chr10:26653933]t       chr10

            let alt_end0 = alt.len() - 1;
            let mut scf_left: Vec<u8>;
            let scf_right: Vec<u8>;
            if alt[alt_end0] == b'[' {
                scf_left = get_ref_seq(faidx,  &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                scf_right = get_ref_seq(faidx, &chrom_2, pos1_2, pos1_2 + SCAFFOLD_FLANK_LEN - 1, false)?;
            } else if alt[alt_end0] == b']' {
                scf_left = get_ref_seq(faidx,  &chrom_1, pos1_1 - SCAFFOLD_FLANK_LEN + 1, pos1_1, false)?;
                scf_right = get_ref_seq(faidx, &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, true)?;
            } else if alt[0] == b']' {
                scf_left = get_ref_seq(faidx,  &chrom_2, pos1_2 - SCAFFOLD_FLANK_LEN + 1, pos1_2, false)?;
                scf_right = get_ref_seq(faidx, &chrom_1, pos1_1, pos1_1 + SCAFFOLD_FLANK_LEN - 1, false)?;
            } else {
                scf_left = get_ref_seq(faidx,  &chrom_2, pos1_2, pos1_2 + SCAFFOLD_FLANK_LEN - 1, true)?;
                scf_right = get_ref_seq(faidx, &chrom_1, pos1_1, pos1_1 + SCAFFOLD_FLANK_LEN - 1, false)?;
            }
            scf_left.extend(scf_right); 
            scfs.push(scf_left);
        },
        _ => {}
    }
    Ok(Some(scfs))
}

/// Fetch a reference sequence given chrom, 1-based start, and 1-based end.
fn get_ref_seq(
    faidx:  &faidx::Reader, 
    chrom:  &str, 
    start1: usize, 
    end1:   usize,
    rc:     bool,
) -> Result<Vec<u8>, Box<dyn Error>> {
    let seq = faidx.fetch_seq(chrom, start1 - 1, end1)?;
    let mut seq = seq.to_ascii_uppercase();
    if rc { 
        seq = seq.iter().rev().map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _    => b'N',
        }).collect(); 
    }
    Ok(seq)
}

// -----------------------------------------------------------------------------
// SEQ extraction
// -----------------------------------------------------------------------------
/// Fill the SEQs for the set of unique QNAMEs found in any SV read into the tool.
fn collect_read_seqs(w: &mut Workflow, tool: &mut DigestionTool) -> Result<(), Box<dyn Error>> {
    w.log.print("extracting read sequences");

    // initialize reading
    let benchmark_bam = w.cfg.get_string(BENCHMARK_BAM).to_string();
    let mut bam = bam::Reader::from_path(benchmark_bam)?;
    bam.set_threads(tool.n_cpu)?;
    let mut aln = bam::Record::new();

    // collect SV-containing read sequences on first encounter in BAM
    while let Some(result) = bam.read(&mut aln) {
        result?;
        match tool.seqs.get_mut(aln.qname()){
            Some(seq) => if seq.is_empty() {
                seq.extend(aln.seq().as_bytes()); // orientation does not matter for mini-consensus
            },
            None => continue
        }
    }
    Ok(())
}

// -----------------------------------------------------------------------------
// site finding in read consensuses
// -----------------------------------------------------------------------------
/// Establish a bcf::Writer updates to the required new tags:
fn get_bcf_writer(tool: &DigestionTool, header: &mut bcf::Header) -> Result<bcf::Writer, Box<dyn Error>>{
    for re in &tool.res {
        let record = format!("##INFO=<ID={},Number=.,Type=String,Description=\"de novo base insertions at junctions\">", SVINSSEQS);
        header.push_record(record.as_bytes());  
        let record = format!("##INFO=<ID={},Number=.,Type=String,Description=\"{} fragment sizes\">", re.size_tag_str, re.name);
        header.push_record(record.as_bytes());  
        let record = format!("##INFO=<ID={},Number=.,Type=String,Description=\"{} nearest site distances\">", re.dist_tag_str, re.name);
        header.push_record(record.as_bytes());  
    }
    let vcf = bcf::Writer::from_stdout(
        header, 
        false, // compress output with BGZF
        bcf::Format::Vcf
    )?;
    Ok(vcf)
}

/// Resolve consensus sequences of SV reads on their scaffolds and perform RE
/// site finding.
fn resolve_consensuses(w: &mut Workflow, tool: &mut DigestionTool, vcf: &mut bcf::Writer) -> Result<(), Box<dyn Error>>{
    w.log.print("resolving consensuses and performing RE site analysis");
    let (mut resolver, mut poa_pool) = Resolver::with_capacity(
        ResolverConfig::default(),
        100, 
        SCAFFOLD_FLANK_LEN * 2 + 1000,
    );

    // process one SV line at a time
    for (sv, scfs) in tool.svs.iter_mut().zip(tool.scfs.iter()) {

        // initialize the SV record for site analysis
        vcf.translate(sv);
        let Some(qnames) = get_vcf_qnames(sv)? else { continue };
        let downsampled_qnames = &qnames[0..qnames.len().min(MAX_CONSENSUS_SEQS)];
        let mut ins_seqs = vec!["?".to_string(); scfs.len()];
        let mut site_metadata: FxHashMap<Vec<u8>, Vec<i32>> = FxHashMap::default();

        // get the fragment sizes for each RE for each junction in the SV line (are two for INV SVs)
        for (scf_i0, scf) in scfs.iter().enumerate() {

            // build the SV allele haplotype consensus sequence
            resolver.set_scaffold_with_aligner(scf);
            for qname in downsampled_qnames {
                let seq = tool.seqs.get(qname).unwrap();
                let _ = resolver.add_seq(seq, |_aln| true);
            }
            let haplotype = resolver.get_consensus(&mut poa_pool);
            let haplotype_str = unsafe{ from_utf8_unchecked(&haplotype) };

            // locate the junction in the consensus and characterize de novo base insertions
            let minimap2 = Minimap2::builder()
                .map_hifi()
                .with_cigar() // required if CS tag is desired
                .with_seq(scf)
                .expect("Failed to build minimap2 index");
            let hap_on_scf = minimap2.map(
                &haplotype, true, false, 
                None, None, None
            )?;
            let mut hap_on_scf_map = vec![None::<usize>; scf.len()];
            for mapping in &hap_on_scf{
                if !mapping.is_primary { break }
                process_cs_tag(&mut hap_on_scf_map, mapping);
            }
            let mut jxn_pos0 = None::<usize>;
            if let (
                Some(left_bkpt0), 
                Some(right_bkpt0)
            ) = (
                hap_on_scf_map[..SCAFFOLD_FLANK_LEN - 1].iter().rev().find_map(|&x| x),
                hap_on_scf_map[scf.len() - SCAFFOLD_FLANK_LEN..].iter().find_map(|&x| x)
            ) {
                jxn_pos0 = Some(left_bkpt0 + (right_bkpt0 - left_bkpt0) / 2);
                ins_seqs[scf_i0] = if right_bkpt0 - left_bkpt0 > 1 {
                    unsafe{ from_utf8_unchecked(&haplotype[left_bkpt0 + 1..right_bkpt0]).to_string() }
                } else {
                    "-".to_string()
                }
            }

            // get fragment metadata for each RE
            for re in &tool.res {
                let (frag_size, min_site_dist) = if let Some(jxn_pos0) = jxn_pos0 {
                    get_site_data(haplotype_str, &re.site, jxn_pos0)
                } else {
                    (-1, -1)
                };
                site_metadata.entry(re.size_tag.clone())
                    .or_insert_with(|| Vec::with_capacity(2))
                    .push(frag_size);
                site_metadata.entry(re.dist_tag.clone())
                    .or_insert_with(|| Vec::with_capacity(2))
                    .push(min_site_dist);
            }
        }

        // add RE fragment metadata INFO tags and write the final VCF record
        let ins_seqs: Vec<_> = ins_seqs.iter().map(|x| x.as_bytes()).collect();
        sv.push_info_string(SVINSSEQS_TAG, &ins_seqs)?;
        for re in &tool.res {
            sv.push_info_integer(&re.size_tag, &site_metadata[&re.size_tag])?;
            sv.push_info_integer(&re.dist_tag, &site_metadata[&re.dist_tag])?;
        }
        vcf.write(sv)?;
    };

    Ok(())
}

/// Process a cs:Z:tag to help fill hap_on_scf_map from one primary mapping.
fn process_cs_tag(
    hap_on_scf_map: &mut Vec<Option<usize>>,
    mapping: &Mapping,
) { 
    let Some(aln) = &mapping.alignment else { return };
    let Some(cs) = &aln.cs else { return };
    let mut hap_pos0 = mapping.query_start as usize;
    let mut scf_pos0 = mapping.target_start as usize;
    let mut chars = cs.chars();
    let mut cs_op = chars.next().unwrap();
    let mut op_val = String::new();
    while let Some(char) = chars.next() {
        if char.is_alphanumeric() {
            op_val.push(char);
        } else {
            handle_cs_op(hap_on_scf_map, cs_op, &op_val, &mut hap_pos0, &mut scf_pos0);
            cs_op = char;
            op_val.clear();
        }
    }
    handle_cs_op(hap_on_scf_map, cs_op, &op_val, &mut hap_pos0, &mut scf_pos0);
}

/// Process one cs:Z:tag operation to add haplotype position to hap_on_scf_map. 
fn handle_cs_op(
    hap_on_scf_map: &mut Vec<Option<usize>>,
    cs_op : char,
    op_val: &str,
    hap_pos0: &mut usize, 
    scf_pos0: &mut usize,
) {
    match cs_op {
        // :[0-9]+   Identical sequence length
        ':' => {
            let op_len = op_val.parse::<usize>().unwrap();
            let mut hap_pos: std::ops::Range<_> = *hap_pos0..*hap_pos0 + op_len;
            hap_on_scf_map[*scf_pos0..*scf_pos0 + op_len].fill_with(|| hap_pos.next());
            *hap_pos0 += op_len;
            *scf_pos0 += op_len;
        },
        // *[acgtn][acgtn]   Substitution: target to query
        '*' => {
            //     S
            // rrrrRrrrr
            // qqqqQqqqq
            //     A
            hap_on_scf_map[*scf_pos0] = Some(*hap_pos0);
            *hap_pos0 += 1;
            *scf_pos0 += 1;
        },
        // +[acgtn]+   Insertion to the target
        '+' => {
            //    *INI     insertions may have heteroduplex bases within homoduplex query run
            // rrrr   Rrrr
            // qqqqQqqqqqq
            //    aA Aa
            let op_len = op_val.len();
            *hap_pos0 += op_len;
        },
        // -[acgtn]+   Deletion from the target
        '-' => {
            //     DDD
            // rrrrRrrrrrr
            // qqqq   Qqqq
            //   aA   Aa
            let op_len = op_val.len();
            hap_on_scf_map[*scf_pos0..*scf_pos0 + op_len].fill(Some(*hap_pos0 - 1));
            *scf_pos0 += op_len;
        },
        _   => panic!("Unexpected operation in cs tag: {}", cs_op),
    }
}

/// Find the RE sites to the left and to the right of the ~center of the junction
/// in the haplotype consensus; return (fragment size, nearest site distance).
fn get_site_data(haplotype_str: &str, site: &str, jxn_pos0: usize) -> (i32, i32) {
    // -------------------------|---------------------------
    //   ######               ######            ######                     
    let half_site_len = site.len() / 2;

    let boundary = (jxn_pos0 + half_site_len).min(haplotype_str.len());
    let left_site_pos = haplotype_str[..boundary]
        .rfind(site)
        .map(|i0| i0 + half_site_len)
        .unwrap_or(0);

    let boundary = jxn_pos0.saturating_sub(half_site_len);
    let mut right_site_pos = haplotype_str[boundary..]
        .find(site)
        .map(|i0| boundary + i0 + half_site_len)
        .unwrap_or(haplotype_str.len());
    if right_site_pos < left_site_pos { right_site_pos = left_site_pos }

    let frag_size = right_site_pos - left_site_pos;
    let left_site_dist  = jxn_pos0.saturating_sub(left_site_pos);
    let right_site_dist = right_site_pos.saturating_sub(jxn_pos0);

    (frag_size as i32, left_site_dist.min(right_site_dist) as i32)
}

// DEL (Deletion)
// DUP (Duplication)
// INV (Inversion)
// BND (Breakend)
// INS (Insertion)

// ##fileformat=VCFv4.2
// ##source=Sniffles2_2.8.1
// ##command="sniffles ...."
// ##fileDate="2026/11/06 18:20:34"
// ##contig=<ID=chr1,length=248387328>
// ...
// ##contig=<ID=chr21,length=45090682>
// ##ALT=<ID=INS,Description="Insertion">
// ##ALT=<ID=DEL,Description="Deletion">
// ##ALT=<ID=DUP,Description="Duplication">
// ##ALT=<ID=INV,Description="Inversion">
// ##ALT=<ID=BND,Description="Breakend; Translocation">
// ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
// ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
// ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
// ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
// ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">
// ##FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">
// ##FILTER=<ID=PASS,Description="All filters passed">
// ##FILTER=<ID=GT,Description="Genotype filter">
// ##FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">
// ##FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">
// ##FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">
// ##FILTER=<ID=COV_MIN,Description="Minimum coverage filter">
// ##FILTER=<ID=COV_MIN_GT,Description="Minimum coverage filter (missing genotype)">
// ##FILTER=<ID=COV_CHANGE_DEL,Description="Coverage change filter for DEL">
// ##FILTER=<ID=COV_CHANGE_DUP,Description="Coverage change filter for DUP">
// ##FILTER=<ID=COV_CHANGE_INS,Description="Coverage change filter for INS">
// ##FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage fractional change filter: upstream-start">
// ##FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage fractional change filter: start-center">
// ##FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage fractional change filter: center-end">
// ##FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage fractional change filter: end-downstream">
// ##FILTER=<ID=COV_VAR,Description="Coverage variance exceeded">
// ##FILTER=<ID=MOSAIC_VAF,Description="Mosaic variant allele fraction filter">
// ##FILTER=<ID=NOT_MOSAIC_VAF,Description="Variant allele fraction filter for non-mosaic">
// ##FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">
// ##FILTER=<ID=STRAND_BND,Description="Strand support filter for BNDs">
// ##FILTER=<ID=STRAND,Description="Strand support filter for germline SVs">
// ##FILTER=<ID=STRAND_MOSAIC,Description="Strand support filter for mosaic SVs">
// ##FILTER=<ID=SVLEN_MIN,Description="SV length filter">
// ##FILTER=<ID=SVLEN_MIN_MOSAIC,Description="SV length filter for mosaic SVs (min)">
// ##FILTER=<ID=SVLEN_MAX_MOSAIC,Description="SV length filter for mosaic SVs (max)">
// ##FILTER=<ID=SINGLE_BREAK,Description="A single break point was detected but not classified as an SV.">
// ##FILTER=<ID=INLINE_SA,Description="INLINE/CIGAR-based SV is mostly supported by SA reads">
// ##FILTER=<ID=MOSAIC_SV_CLOSE_EDGE,Description="For mosaic SVs, the location is close to the end of the read (either end)">
// ##FILTER=<ID=GT_FAILED,Description="Sniffles was unable to genotype this call.">
// ##FILTER=<ID=UNRESOLVED,Description="Call is within an unresolved region.">
// ##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">
// ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">
// ##INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Structural variation classified as putative mosaic">
// ##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">
// ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">
// ##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">
// ##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">
// ##INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only), or reads for a DUP SV reported as INS by the mapper.">
// ##INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of reads supporting a DEL SV through supplementary alignments (split events)">
// ##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">
// ##INFO=<ID=SUPPORT_UNSCALED,Number=1,Type=Float,Description="Base support for long insertions, before adjustment">
// ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">
// ##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">
// ##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">
// ##INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">
// ##INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">
// ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">
// ##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">
// ##INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">
// ##INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">
// ##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">
// ##INFO=<ID=COVERAGE_VAR,Number=1,Type=Float,Description="Variance of coverage across large events">
// ##INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">
// ##INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">
// ##INFO=<ID=LASM,Number=0,Type=Flag,Description="Local assembly used to detect the structural variant">
// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE

// chr1    2707    Sniffles2.DEL.53DS0     TCTgca...ccg      T       60      PASS    PRECISE;SVTYPE=DEL;SVLEN=-259;END=2966;SUPPORT=4;RNAMES=m64049_210425_091619/25821819/ccs,m64018_210503_204922/153747606/ccs,m64049_210424_024840/112986691/ccs,m64018_210505_030429/14092581/ccs;COVERAGE=5,5,16,21,22;STRAND=+-;NM=0.027;PHASE=0,NULL,4,4,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;SUPPORT_SA=0;VAF=0.286        GT:GQ:DR:DV:PS  0/1:12:10:4:.

// chr1    3009031 Sniffles2.DUP.B92S0     G       <DUP>   60      PASS    IMPRECISE;SVTYPE=DUP;SVLEN=1394;END=3010425;SUPPORT=20;RNAMES=m64109_200807_075817/122487993/ccs,m64049_210425_091619/157156687/ccs,m64018_210505_030429/147261366/ccs,m64049_210424_024840/14484519/ccs,m64049_210424_024840/28770683/ccs,m64049_210424_024840/145623701/ccs,m64018_210505_030429/49219233/ccs;COVERAGE=71,87,85,85,70;STRAND=+-;NM=0.008;PHASE=0,NULL,7,7,PASS,FAIL;STDEV_LEN=17.964;STDEV_POS=14.789;SUPPORT_INLINE=13;VAF=0.286 GT:GQ:DR:DV:PS  0/1:60:50:20:.

// chr1    16396440        Sniffles2.INV.BBBS0     C       <INV>   59      PASS    PRECISE;SVTYPE=INV;SVLEN=128772785;END=145169225;SUPPORT=57;RNAMES=m64109_200805_204709/125502700/ccs,m64018_210503_204922/24182989/ccs,...;COVERAGE=75,137,66,76,138;STRAND=+;NM=0.019;PHASE=0,NULL,57,57,PASS,FAIL;STDEV_LEN=0.000;STDEV_POS=0.000;VAF=0.425    GT:GQ:DR:DV:PS    0/1:60:77:57:.

// chr1    16396441        Sniffles2.BND.BDCS0     C       [chr1:145169227[C       59      PASS    PRECISE;SVTYPE=BND;SUPPORT=60;RNAMES=m64109_200805_204709/125502700/ccs,m64018_210503_204922/24182989/ccs,...;COVERAGE=75,75,138,140,143;STRAND=+-;NM=606.417;CHR2=chr1;PHASE=0,NULL,60,60,PASS,FAIL;STDEV_POS=0.000;VAF=0.508       GT:GQ:DR:DV:PS  0/1:60:58:60:.

// chr1    2566    Sniffles2.INS.0S0       c       cCCTCACCCTCACCCTCACCCTCACCCTGACCCTGACCCTGACCCTCACCCTGACCCTCACCCTCACCCTCA        60      PASS    IMPRECISE;SVTYPE=INS;SVLEN=71;END=2566;SUPPORT=3;RNAMES=m64018_210505_030429/14092581/ccs,m64018_210503_204922/153747606/ccs,...;COVERAGE=5,5,5,5,22;STRAND=+-;NM=0.026;PHASE=0,NULL,3,3,PASS,FAIL;STDEV_LEN=5.033;STDEV_POS=54.271;SUPPORT_LONG=0;SUPPORT_UNSCALED=3;VAF=0.600      GT:GQ:DR:DV:PS  0/1:11:2:3:.
