//! Support for analyzing PacBio kinetics tags in BAM records
//! toward assessing base damage when strand base values differ.

// dependencies
use crossbeam::channel::Sender;
use crate::tools::basecall_pacbio::{BufferedStrand, KineticsInstance};

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
pub struct StrandMetadata {
    pub base_context:    [u8; 3], // 3-base context around the query base
    pub is_heteroduplex: bool,    // only is_heteroduplex is set on initial basecalling; is_ref fields set after alignment
    pub one_is_ref:      bool,    // whether one of the strands is the reference base
    pub this_is_ref:     bool,    // whether this strand is the reference base
}

/// FrameCounts structure for storing the pulse values 
/// surrounding a query base position. 
pub struct FrameCounts {
    pub ip_before:   u16,  // inter-pulse duration before base incorporation
    pub pw:          u16,  // pulse width at base incorporation
    pub ip_after:    u16,  // inter-pulse duration after base incorporation
}

/// Extract the FrameCounts for a base at the given offset in the strand.
/// Send a KineticsInstance to the kinetics collector via tx_kinetics.
pub fn push(
    strand: &BufferedStrand,
    offset: usize,
    is_heteroduplex: bool,
    tx_kinetics: &Sender<KineticsInstance>,
) -> (String, Vec<u16>) {

    // ignore identity bases at the ends of the read 
    if offset == 0 || offset + 1 >= strand.seq.len() {
        return ("NNN".to_string(), vec![0, 0, 0]); // dummy value, never expeced to be used in caller
    }

    // extract the 3-base context around the target base
    let base_context = &strand.seq[offset - 1..=offset + 1];

    // extract and decode kinetics values
    // do nothing if kinetics data is missing
    let (ip_before, pw, ip_after) = if let (Some(ip), Some(pw)) = (&strand.ip, &strand.pw) {
        (
            CodecV1::decode(ip[offset]),     // inter-pulse duration leading into base incorporation
            CodecV1::decode(pw[offset]),     // pulse width at base incorporation
            CodecV1::decode(ip[offset + 1]), // inter-pulse duration after base incorporation
        )
    } else {
        return (base_context.to_string(), vec![0, 0, 0]);
    };

    // assemble the StrandMetadata key, here based on base_context and is_heteroduplex only
    let strand_metadata = StrandMetadata {
        base_context: base_context.as_bytes().try_into().unwrap(),
        is_heteroduplex,
        one_is_ref:   false,
        this_is_ref:  false,
    };

    // assemble the decoded FrameCounts
    let frame_counts = FrameCounts {
        ip_before,
        pw,
        ip_after,
    };

    // send the kinetics instance to the collector for writing
    tx_kinetics.send(KineticsInstance {
        strand_metadata,
        frame_counts,
    }).unwrap();

    // return the 3-base context string for dd:Z: tag
    (base_context.to_string(), vec![ip_before, pw, ip_after])
}
