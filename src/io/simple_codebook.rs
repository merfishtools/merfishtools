use std::path::Path;

use itertools::Itertools;
use serde::de::Error;
use serde::de::Unexpected;
use serde::{Deserialize, Deserializer, Serialize};

use crate::simulation::binary;

/// Deserialize bool from String with custom value mapping
fn bool_from_string<'de, D>(deserializer: D) -> Result<bool, D::Error>
    where
        D: Deserializer<'de>,
{
    match String::deserialize(deserializer)?.as_ref() {
        "1" => Ok(true),
        "0" => Ok(false),
        other => Err(Error::invalid_value(Unexpected::Str(other), &"1 or 0")),
    }
}

/// A codebook record.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Record {
    #[serde(rename = "feat")]
    name: String,
    #[serde(with = "binary")]
    codeword: u16,
    #[serde(deserialize_with = "bool_from_string")]
    expressed: bool,
}

impl Record {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn codeword(&self) -> u16 {
        self.codeword
    }

    pub fn expressed(&self) -> bool {
        self.expressed
    }
}

pub struct SimpleCodebook {
    records: Vec<Record>,
}

impl SimpleCodebook {
    pub fn from_file<P: AsRef<Path>>(path: P) -> csv::Result<Self> {
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
        let records: Vec<Record> = rdr
            .deserialize()
            .map(|r| r.expect("Failed reading codebook record."))
            .collect_vec();
        Ok(SimpleCodebook { records })
    }
    pub fn records(&self) -> &Vec<Record> {
        &self.records
    }

    pub fn contains(&self, name: &str) -> bool {
        self.records.iter().any(|r| r.name() == name)
    }

    pub fn num_bits(&self) -> u16 {
        let max_cw = self.records.iter().map(|r| r.codeword).max();
        if let Some(bits) = max_cw {
            (bits as f64).log2().ceil() as u16
        } else {
            16
        }
    }

    pub fn num_bits_set(&self) -> u16 {
        let max_bits_set = self.records.iter().map(|r| r.codeword.count_ones()).max();
        if let Some(bits) = max_bits_set {
            bits as u16
        } else {
            4
        }
    }
}
