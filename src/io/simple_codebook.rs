use std::path::Path;

use itertools::Itertools;
use serde::{Deserialize, Deserializer, Serialize};
use serde::de::Error;
use serde::de::Unexpected;

use crate::simulation::binary;

/// Deserialize bool from String with custom value mapping
fn bool_from_string<'de, D>(deserializer: D) -> Result<bool, D::Error>
    where
        D: Deserializer<'de>,
{
    match String::deserialize(deserializer)?.as_ref() {
        "1" => Ok(true),
        "0" => Ok(false),
        other => Err(Error::invalid_value(
            Unexpected::Str(other),
            &"1 or 0",
        )),
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

    pub fn codeword(&self) -> u16 { self.codeword }

    pub fn expressed(&self) -> bool {
        self.expressed
    }
}

pub struct SimpleCodebook {
    records: Vec<Record>
}

impl SimpleCodebook {
    pub fn from_file<P: AsRef<Path>>(path: P) -> csv::Result<Self> {
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
        let records: Vec<Record> = rdr.deserialize().map(|r| r.expect("Failed reading codebook record.")).collect_vec();
        Ok(SimpleCodebook {
            records
        })
    }
    pub fn records(&self) -> &Vec<Record> {
        &self.records
    }

    pub fn contains(&self, name: &str) -> bool {
        self.records.iter().any(|r| r.name() == name)
    }
}