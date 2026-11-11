//! VCF handling during benchmarking.

// imports
use std::error::Error;
use std::str::from_utf8_unchecked;
use rust_htslib::bcf;
use crate::tools::compare_sv_calls::ComparisonTool;
use super::*;

/// One of four possible orientations for BND SV calls.
pub enum BndType {
    /// t[p[  piece extending to the right of p is joined after t
    AfterRight,     
    /// t]p]  reverse comp piece extending left of p is joined after t
    AfterRcLeft,    
    /// ]p]t  piece extending to the left of p is joined before t
    BeforeLeft,     
    /// [p[t  reverse comp piece extending right of p is joined before t
    BeforeRcRight,  
}
impl BndType {
    pub fn from_alt(alt: &[u8]) -> Result<BndType, Box<dyn Error>> {
        let max_i0 = alt.len() - 1;
        if alt[max_i0] == b'[' { return Ok(BndType::AfterRight) }
        if alt[max_i0] == b']' { return Ok(BndType::AfterRcLeft) }
        if alt[0]      == b']' { return Ok(BndType::BeforeLeft) }
        if alt[0]      == b'[' { return Ok(BndType::BeforeRcRight) }
        Err("Failed to parse BndType".into())
    }
}
/// Structure carrying unpacker ALT encoding.
pub struct BndAlt {
    chrom: String,
    pos1:  u32,
    _bnd_type: BndType
}

// -----------------------------------------------------------------------------
// VCF path extraction
// -----------------------------------------------------------------------------
/// Get VCF ALT as a (BND packed) String.as_bytes().
pub fn get_vcf_alt(jxn: &bcf::Record) -> Result<&[u8], Box<dyn Error>> {
    let alleles = jxn.alleles();
    let alt  = *alleles.get(1).ok_or("Missing ALT")?;
    Ok(alt)
}

/// Get VCF left breakpoint node.
pub fn get_vcf_left_node(jxn: &bcf::Record) -> Result<(String, u32), Box<dyn Error>> {
    let chrom_idx_1 = jxn.rid().ok_or("Missing RID")?;
    let chrom_bytes_1 = jxn.header().rid2name(chrom_idx_1)?;
    let chrom_1 = unsafe { from_utf8_unchecked(chrom_bytes_1) };
    let pos1_1 = (jxn.pos() + 1) as u32;
    Ok((chrom_1.to_string(), pos1_1))
}

/// Get VCF right breakpoint node.
pub fn get_vcf_right_node(
    jxn: &bcf::Record, 
    chrom_1: &str, 
    alt:     &[u8], 
    sv_type: &[u8]
) -> Result<(String, u32), Box<dyn Error>> {
    let chrom_2 = match jxn.info(b"CHR2").string()? {
        Some(x) => unsafe { // followed by Sniffles2 BND
            from_utf8_unchecked(x[0]).to_string() // may or may not be chrom_1
        },
        None => match sv_type {
            b"BND" => unpack_alt(alt)?.chrom, // followed by CuteSV BND
            _ => chrom_1.to_string(),
        }
    };
    let pos1_2 = match jxn.info(b"END").integer()? {
        Some(v) => (v[0] as u32).to_owned(),
        None => { // BND path; callers do not use END for BND
            unpack_alt(alt)?.pos1
        },
    };
    Ok((chrom_2.to_string(), pos1_2))
}
fn unpack_alt(alt: &[u8]) -> Result<BndAlt, Box<dyn Error>> {
    let _bnd_type = BndType::from_alt(alt)?;
    let alt =  unsafe { from_utf8_unchecked(alt) };
    let first_bracket = alt.find('[').or_else(|| alt.find(']'))
        .ok_or("Missing first bracket")?;
    let rest = &alt[first_bracket + 1..];
    let second_bracket = rest.find('[').or_else(|| rest.find(']'))
        .ok_or("Missing second bracket")?;
    let bnd_target = &rest[..second_bracket]; // e.g., "chr1:145169227"
    let (chrom, pos1) = bnd_target.split_once(':')
        .ok_or("Bad chr:pos in ALT")?;
    Ok(BndAlt {
        chrom: chrom.to_string(),
        pos1:  pos1.parse::<u32>()?,
        _bnd_type
    })
}

// -----------------------------------------------------------------------------
// VCF INFO extraction
// -----------------------------------------------------------------------------
/// Get INFO tag SVTYPE as Option<Vec<u8>>.
pub fn get_vcf_sv_type(jxn: &bcf::Record) -> Result<Option<Vec<u8>>, Box<dyn Error>> {
    match jxn.info(b"SVTYPE").string()? {
        Some(x) => Ok(Some(x[0].to_vec())),
        None => Ok(None)
    }
}
/// Get INFO tag SVTYPE as String; panics if SVTYPE not found.
pub fn get_vcf_sv_type_str(jxn: &bcf::Record) -> String {
    let Ok(Some(sv_type)) = get_vcf_sv_type(jxn) else {
        panic!("Missing SVTYPE on junction");
    };
    unsafe{ from_utf8_unchecked(&sv_type).to_string() }
}

/// Get INFO tag STRAND (CuteSV) or "**" (Sniffles, its STRAND field is for reads). 
pub fn get_vcf_strand_bytes(jxn: &bcf::Record, source: &SvSource) -> Result<Vec<u8>, Box<dyn Error>> {
    match source {
        SvSource::Sniffles => Ok(b"**".to_vec()),
        _ => match jxn.info(b"STRAND").string()? {
            Some(x) => Ok(x[0].to_vec()),
            None => Ok(b"**".to_vec()),
        },
    }
}
/// Get INFO tag STRAND (CuteSV) or "**" (Sniffles, its STRAND field is for reads). 
pub fn get_vcf_strand_str(jxn: &bcf::Record, source: &SvSource) -> Result<String, Box<dyn Error>> {
    let strand = get_vcf_strand_bytes(jxn, source)?;
    unsafe{ Ok(from_utf8_unchecked(&strand).to_string()) }
}

/// Get the VCF SV size from chrom and pos, using 999_999_999 for translocations.
pub fn get_vcf_sv_size(chrom_1: &str, pos1_1: u32, chrom_2: &str, pos1_2: u32) -> u32{
    if chrom_1 == chrom_2 {
        (pos1_2 as i32 - pos1_1 as i32).unsigned_abs()
    } else {
        999_999_999
    }
}

/// Get RE fragment size added by digest_sv_alleles.
pub fn get_frag_size(jxn: &bcf::Record, tool: &ComparisonTool) -> Result<i32, Box<dyn Error>> {
    match jxn.info(&tool.re.size_tag).integer() {
        Ok(Some(v)) => Ok(v[0].to_owned()),
        _ => Ok(0) // the path taken when RE is NA
    }
}

/// Get nearest site distance added by digest_sv_alleles.
pub fn get_stem_length(jxn: &bcf::Record, tool: &ComparisonTool) -> Result<u32, Box<dyn Error>> {
    match jxn.info(&tool.re.dist_tag).integer() {
         Ok(Some(v)) => Ok((v[0] as u32).to_owned()),
        _ => Ok(0) // the path taken when RE is NA
    }
}

/// Get INFO tag SVINSSEQ (as added by digest_sv_alleles).
pub fn get_vcf_ins_seq(jxn: &bcf::Record) -> Result<String, Box<dyn Error>> {
    match jxn.info(SVINSSEQ_TAG).string()? {
        Some(x) => Ok(
            unsafe{ from_utf8_unchecked(x[0]).to_string() }
        ),
        None => Ok("?".to_string()),
    }
}

/// Get INFO tag SUPPORT (Sniffles) or RE (CuteSV).
pub fn get_vcf_n_observed(jxn: &bcf::Record, source: &SvSource) -> Result<u16, Box<dyn Error>> {
    let info = match source {
        SvSource::Sniffles => jxn.info(b"SUPPORT"),
        _ => jxn.info(b"RE"),
    };
    match info.integer()? {
        Some(x) => Ok((x[0] as u16).to_owned()),
        None => Ok(0),
    }
}

/// Get INFO tag MOSAIC (Sniffles) or `false` (CuteSV).
pub fn get_vcf_mosaic(jxn: &bcf::Record, source: &SvSource) -> Result<bool, Box<dyn Error>> {
    match source {
        SvSource::CuteSV => Ok(false),
        _ => match jxn.info(b"MOSAIC").flag()? {
            x=> Ok(x)
        },
    }
}

/// Get INFO tag VAF (Sniffles) or AF (CuteSV).
pub fn get_vcf_vaf(jxn: &bcf::Record, source: &SvSource) -> Result<f32, Box<dyn Error>> {
    let info = match source {
        SvSource::Sniffles => jxn.info(b"VAF"),
        _ => jxn.info(b"AF"),
    };
    match info.float()? {
        Some(x) => Ok((x[0] as f32).to_owned()),
        None => Ok(0.0),
    }
}

/// Extract the QNAMES from INFO/RNAMES (which should be called QNAMES).
pub fn get_vcf_qnames(jxn: &bcf::Record) -> Result<Option<Vec<Vec<u8>>>, Box<dyn Error>> {
    match jxn.info(b"RNAMES").string()? {
        Some(x) => {
            Ok(Some(x.iter().map(|&qname| qname.to_vec()).collect()))
        },
        None => Ok(None),
    }
}
