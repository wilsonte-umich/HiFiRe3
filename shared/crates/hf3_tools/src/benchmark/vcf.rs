//! VCF handling during benchmarking.

// imports
use std::error::Error;
use std::str::from_utf8_unchecked;
use rust_htslib::bcf;
use super::*;

// -----------------------------------------------------------------------------
// VCF path extraction
// -----------------------------------------------------------------------------
/// Get VCF ALT.
pub fn get_vcf_alt(sv: &bcf::Record) -> Result<&[u8], Box<dyn Error>> {
    let alleles = sv.alleles();
    let alt  = *alleles.get(1).ok_or("Missing ALT")?;
    Ok(alt)
}

/// Get VCF left breakpoint node.
pub fn get_vcf_left_node(sv: &bcf::Record) -> Result<(String, usize), Box<dyn Error>> {
    let chrom_idx_1 = sv.rid().ok_or("Missing RID")?;
    let chrom_bytes_1 = sv.header().rid2name(chrom_idx_1)?;
    let chrom_1 = unsafe { from_utf8_unchecked(chrom_bytes_1) };
    let pos1_1 = (sv.pos() + 1) as usize;
    Ok((chrom_1.to_string(), pos1_1))
}

/// Get VCF right breakpoint node.
pub fn get_vcf_right_node(sv: &bcf::Record, alt: &[u8], chrom_1: &str) -> Result<(String, usize), Box<dyn Error>> {
    let chrom_2 = match sv.info(b"CHR2").string()? {
        Some(x) => unsafe { 
            from_utf8_unchecked(x[0]).to_string() // may or may not be chrom_1
        },
        None => chrom_1.to_string(), // all but BND
    };
    let pos1_2 = match sv.info(b"END").integer()? {
        Some(v) => (v[0] as usize).to_owned(),
        None => { // BND path; Sniffles2 does not use END for BND
            let alt =  unsafe { from_utf8_unchecked(alt) };
            let first_bracket = alt.find('[').or_else(|| alt.find(']'))
                .ok_or("Missing first bracket")?;
            let rest = &alt[first_bracket + 1..];
            let second_bracket = rest.find('[').or_else(|| rest.find(']'))
                .ok_or("Missing second bracket")?;
            let bnd_target = &rest[..second_bracket]; // e.g., "chr1:145169227"
            let (_chrom, pos1) = bnd_target.split_once(':')
                .ok_or("Bad chr:pos in ALT")?;
            pos1.parse::<usize>()?
        },
    };
    Ok((chrom_2.to_string(), pos1_2))
}

// -----------------------------------------------------------------------------
// VCF INFO extraction
// -----------------------------------------------------------------------------
/// Get INFO tag SVTYPE.
pub fn get_vcf_sv_type(sv: &bcf::Record) -> Result<Option<Vec<u8>>, Box<dyn Error>> {
    match sv.info(b"SVTYPE").string()? {
        Some(x) => Ok(Some(x[0].to_vec())),
        None => Ok(None)
    }
}

/// Get INFO tag SVLEN.
pub fn get_vcf_sv_len(sv: &bcf::Record) -> Result<Option<u32>, Box<dyn Error>> {
    match sv.info(b"SVLEN").integer()? {
        Some(x) => Ok(Some(x[0].unsigned_abs())),
        None => Ok(None),
    }
}

/// Get INFO tag SVLEN.
pub fn get_vcf_sv_len_with_default(sv: &bcf::Record) -> Result<u32, Box<dyn Error>> {
    match sv.info(b"SVLEN").integer()? {
        Some(x) => Ok(x[0].unsigned_abs()),
        None => Ok(0),
    }
}

/// Extract the QNAMES from INFO/RNAMES (which should be called QNAMES).
pub fn get_vcf_qnames(sv: &bcf::Record) -> Result<Option<Vec<Vec<u8>>>, Box<dyn Error>> {
    match sv.info(b"RNAMES").string()? {
        Some(x) => {
            Ok(Some(x.iter().map(|&qname| qname.to_vec()).collect()))
        },
        None => Ok(None),
    }
}

/// Get INFO tag SUPPORT_SA.
pub fn get_vcf_n_supplemental(sv: &bcf::Record) -> Result<u16, Box<dyn Error>> {
    match sv.info(b"SUPPORT_SA").integer()? {
        Some(x) => Ok((x[0] as u16).to_owned()),
        None => Ok(0),
    }
}

/// Get INFO tag SVINSSEQS (as added by digest_sv_allele).
pub fn get_vcf_ins_seqs(sv: &bcf::Record, i: usize) -> Result<String, Box<dyn Error>> {
    match sv.info(SVINSSEQS_TAG).string()? {
        Some(x) => Ok(
            unsafe{ from_utf8_unchecked(x[i]).to_string() }
        ),
        None => Ok("MISING".to_string()),
    }
}


