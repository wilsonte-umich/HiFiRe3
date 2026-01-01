//! Support for analyzing PacBio kinetics tags in BAM records
//! toward assessing base damage when strand base values differ.

// dependencies
use std::fs::File;
use std::io::Write;
use std::str::from_utf8_unchecked;
use std::collections::HashMap;

// constants
const MAX_CONTEXT_VALUES: usize = 1_000_000; // pre-allocate up to one million entries per base context

/// CodecV1 structure for encoding/decoding for PacBio kinetics values.
struct CodecV1;
impl CodecV1 {
    /// Decode a PacBio CodecV1 byte (u8) into a frame count (u16).
    pub fn decode(b: u8) -> u16 {
        match b {
              0..= 63 =>        b as u16,
             64..=127 =>  64 + (b as u16 -  64) * 2,
            128..=191 => 192 + (b as u16 - 128) * 4,
            192..=255 => 448 + (b as u16 - 192) * 8,
        }
    }
}
// The lossy encoding for IPD and pulsewidth values into the available 
// 256 codepoints is as follows (codec v1):
// Frames	        Encoding
// 0 .. 63	         0, 1, .. 63
// 64, 66, .. 190	 64, 65, .. 127  // 65 is missing since the division and rounding cannot produce it...
// 192, 196 .. 444	 128, 129 .. 191
// 448, 456, .. 952  192, 193 .. 255
// In other words, we use the first 64 codepoints to encode frame counts at single frame resolution, the next 64 to encode the frame counts at two-frame resolution, and so on. Durations exceeding 952 frames are capped at 952. Durations not enumerated in “Frames” above are rounded to the nearest enumerated duration then encoded. For example, a duration of 194 frames would round to 196 and then be encoded as codepoint 129

/// StrandMetadata structure for statifying strand kinetics values
/// by duplex and reference base status at a central query base.
#[derive(Hash, Eq, PartialEq)]
struct StrandMetadata {
    base_context:    [u8; 3], // 3-base context around the query base
    is_heteroduplex: bool,    // only is_heteroduplex is set on initial basecalling; is_ref fields set after alignment
    one_is_ref:      bool,    // whether one of the strands is the reference base
    this_is_ref:     bool,    // whether this strand is the reference base
}

/// FrameCounts structure for storing the pulse values 
/// surrounding a query base position. 
struct FrameCounts {
    ip_before:   u16,  // inter-pulse duration before base incorporation
    pw:          u16,  // pulse width at base incorporation
    ip_after:    u16,  // inter-pulse duration after base incorporation
}

/// Kinetics structure for caching PacBio kinetics values
/// stratified by 3-base context and strand metadata.
pub (super) struct Kinetics {
    map: HashMap<StrandMetadata, Vec<FrameCounts>>, // multiple kinetics values trios per 3-base context and strand type
}
impl Kinetics {
    /// Create a new, empty Kinetics structure.
    pub fn new() -> Self {
        Kinetics {
            map: HashMap::new(),
        }
    }
    /// Cache a set of kinetics values for a 3-base context key,
    /// up to a maximum of MAX_CONTEXT_VALUES per key.
    pub fn push(
        &mut self, 
        seq:    &str, // in the same orientation as the original strand and thus as ip and pw
        ip:     &Option<Vec<u8>>, 
        pw:     &Option<Vec<u8>>, 
        offset: usize, // criteria in the calling function ensure offset is valid
        is_heteroduplex: bool
    ) -> String {

        // extract the 3-base context around the target base
        let base_context = &seq[offset - 1..=offset + 1];

        // extract and decode kinetics values
        // do nothing if kinetics data is missing
        let (ip_before, pw, ip_after) = if let (Some(ip), Some(pw)) = (ip, pw) {
            (
                ip[offset],     // inter-pulse duration leading into base incorporation
                pw[offset],     // pulse width at base incorporation
                ip[offset + 1], // inter-pulse duration after base incorporation
            )
        } else {
            return base_context.to_string();
        };

        // assemble the StrandMetadata key, here based on base_context and is_heteroduplex only
        let strand_metadata = StrandMetadata {
            base_context: base_context.as_bytes().try_into().unwrap(),
            is_heteroduplex,
            one_is_ref:   false,
            this_is_ref:  false,
        };

        // insert into the vector for this context key up to the maximum allowed entries
        let entry = self.map
            .entry(strand_metadata)
            .or_insert(Vec::with_capacity(MAX_CONTEXT_VALUES));
        if entry.len() < MAX_CONTEXT_VALUES {
            let frame_counts = FrameCounts {
                ip_before: CodecV1::decode(ip_before),
                pw:        CodecV1::decode(pw),
                ip_after:  CodecV1::decode(ip_after),
            };
            entry.push(frame_counts);
        }
        base_context.to_string()
    }
    /// Write the cached kinetics values to a PACBIO_BASECALL_KINETICS file.
    pub fn write_kinetics_file(
        &self,
        file_path: &str,
    ) -> Result<(), Box<dyn std::error::Error>> { 
        let mut file = File::create(file_path)?;
        writeln!(
            file,
            "{},{},{},{},{},{},{}",
            "base_context",
            "is_heteroduplex",
            "one_is_ref",
            "this_is_ref",
            "ipd_before",
            "pulse_width",
            "ipd_after"
        )?;
        for (metadata, frame_counts_vec) in &self.map {
            let base_context = unsafe { from_utf8_unchecked(&metadata.base_context) };
            for frame_counts in frame_counts_vec {
                writeln!(
                    file,
                    "{},{},{},{},{},{},{}",
                    base_context,
                    metadata.is_heteroduplex as u8,
                    metadata.one_is_ref as u8,
                    metadata.this_is_ref as u8,
                    frame_counts.ip_before,
                    frame_counts.pw,
                    frame_counts.ip_after,
                )?;
            }
        }
        Ok(())
    }
}
