//! Module for parsing and handling SAM alignment records.

// dependencies
use std::str::FromStr;
use std::fmt;
use serde::{de::{self, Deserializer, Visitor, SeqAccess}, Deserialize, Serialize};
use regex::Regex;

// modules
pub mod flag;

// structures
#[derive(Serialize, Deserialize)]
pub struct SamRecord { // a single SAM alignment record
    pub qname: String,
    pub flag:  u16,
    pub rname: String,
    pub pos1:  u32, 
    pub mapq:  u8,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen:  i32,
    pub seq:   String, // when streamed from STDIN, cannot used borrowed fields, must be DeserializeOwned
    pub qual:  String,
    #[serde(deserialize_with = "SamRecord::deserialize_tags")]
    pub tags:  Vec<String>,
}
impl SamRecord {

    pub fn deserialize_tags<'de, D>(deserializer: D) -> Result<Vec<String>, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct TagsVisitor;
        impl<'de> Visitor<'de> for TagsVisitor {
            type Value = Vec<String>;
            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("row with at least 11 fields")
            }
            fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                // // Skip first 11 fields (indices 0-10: qname through qual)
                // for _ in 0..11 {
                //     let _: serde::de::IgnoredAny = seq.next_element()?.ok_or_else(|| 
                //         de::Error::custom("fewer than 11 fields")
                //     )?;
                // }
                // // Collect remaining fields
                // let mut tags = Vec::new();
                // while let Some(result) = seq.next_element::<String>()? {
                //     tags.push(result);
                // }
                // Ok(tags)
                let mut tags = Vec::new();
                while let Some(tag) = seq.next_element::<String>()? {
                    tags.push(tag);
                }
                Ok(tags)
            }
        }
        deserializer.deserialize_seq(TagsVisitor)
    }

    // check a SAM flag against one or more flag values
    pub fn check_flag(&self, mask: &u16) -> bool {
        (self.flag & mask) != 0
    }

    // extract a parsed tag value from the SAM record using regex
    // returns None if the tag is not present
    pub fn get_tag<T>(&self, tag: &str) -> Result<Option<T>, Box<dyn std::error::Error>> 
    where 
        T: FromStr, 
        <T as FromStr>::Err: fmt::Debug,
    {
        let pattern = format!(r"{}:[^:]+:([^\t]+)", tag);
        let re = Regex::new(&pattern).unwrap();
        for tag in &self.tags {
            if let Some(caps) = re.captures(tag) {
                let value_str = &caps[1];
                match value_str.parse::<T>() {
                    Ok(val) => return Ok(Some(val)),
                    Err(e) => return Err(format!("Failed to parse tag value: {:?}", e).into()),
                }
            }
        }
        Ok(None)
    }
}
