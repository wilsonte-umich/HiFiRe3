
/// BufferedStrand holds the data for a single PacBio strand while it waits
/// for its matching strand to be encountered for merging.
pub struct BufferedStrand {
    pub usable:   bool,
    pub ec:       f32,
    pub seq:      String,
    pub qual:     Vec<u8>,
    pub ip:       Vec<u8>, // inter-pulse durations
    pub pw:       Vec<u8>, // pulse widths
}
