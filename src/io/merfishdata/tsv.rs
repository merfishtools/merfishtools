use core::borrow::BorrowMut;
use std::fs;
use std::io;
use std::iter::Map;
use std::path::Path;

use csv;
use csv::DeserializeRecordsIter;
use rand::distributions::WeightedIndex;
use rand::prelude::StdRng;
use rand::Rng;

use crate::io::codebook::Codebook;
use crate::io::counts::{CommonRecord, FromRecord, into_u16};
use crate::io::merfishdata::{MerfishRecord, Readout};

/// A 2D position in the microscope.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Position {
    pub x: f32,
    pub y: f32,
}

/// A MERFISH raw data record.
/// // "Cell_ID	Gene_Name	Hamming_Distance	Cell_Position_X	Cell_Position_Y	RNA_Position_X	RNA_Position_Y
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TsvRecord {
    #[serde(rename = "cell")]
    pub cell_id: String,
    #[serde(rename = "feat")]
    pub feature: String,
    #[serde(rename = "dist")]
    pub hamming_dist: u8,
    #[serde(flatten, with = "prefix_cell")]
    pub cell_position: Position,
    #[serde(flatten, with = "prefix_none")]
    pub rna_position: Position,
    #[serde(skip)]
    pub codeword: Option<u16>,
    #[serde(skip)]
    pub readout: u16,
}

with_prefix!(prefix_cell "cell_");
with_prefix!(prefix_none "");

impl MerfishRecord for TsvRecord {
    fn cell_id(&self) -> u32 {
        self.cell_id.parse().expect("Failed parsing cell_id")
    }

    fn cell_name(&self) -> String {
        self.cell_id.clone()
    }

    fn cell_pos(&self) -> (f32, f32) {
        (self.cell_position.x, self.cell_position.y)
    }

    fn feature_id(&self) -> u16 {
        self.feature.parse().expect("Failed parsing feature_id")
    }

    fn feature_name(&self) -> String {
        self.feature.clone()
    }

    fn hamming_dist(&self) -> u8 {
        self.hamming_dist
    }

    fn error_mask(&self) -> u16 {
        if let Some(codeword) = self.codeword {
            codeword ^ self.readout
        } else {
            unimplemented!()
        }
    }

    fn codeword(&self) -> Option<u16> {
        self.codeword
    }

    fn readout(&self) -> u16 {
        self.readout
    }

    fn readout_bitvec(&self) -> Readout {
        unimplemented!()
    }

    fn count(&self) -> usize {
        1
    }

    fn is_exact(&self) -> bool {
        self.hamming_dist == 0
    }
}

/// A reader for MERFISH raw data.
pub struct TsvReader<R: io::Read> {
    inner: csv::Reader<R>,
}

pub(crate) fn modify_record(
    r: &TsvRecord,
    codebook: &Codebook,
    num_bits: usize,
    p0: &[f32],
    p1: &[f32],
    rng: &mut StdRng,
) -> TsvRecord {
    let unif = rand::distributions::Uniform::new(0f32, 1.);
    let mut r = r.clone();
    r.codeword = Some(into_u16(
        codebook
            .record(codebook.get_id(&r.feature_name()))
            .codeword(),
    ));
    let codeword = r.codeword.unwrap();
    r.readout = codeword;
    if r.hamming_dist() == 1 {
        let weights: Vec<_> = (0..num_bits)
            .map(|i| {
                let bit = (codeword >> i) & 1;
                match bit {
                    0 => p0[i],
                    1 => p1[i],
                    _ => panic!("bug: bit can only be 0 or 1"),
                }
            })
            .collect();
        // possibly introduce a one-bit error according to weights.
        if rng.sample(unif) > weights.iter().map(|p| 1. - p).product() {
            let wc = WeightedIndex::new(weights).unwrap();
            let error_bit = rng.sample(&wc) as u8;
            r.readout = codeword ^ (1 << error_bit as u16);
        }
    }
    r
}

impl TsvReader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(TsvReader::new)
    }
}

impl<R: io::Read> TsvReader<R> {
    pub fn new(rdr: R) -> Self {
        TsvReader {
            inner: csv::ReaderBuilder::new().delimiter(b'\t').from_reader(rdr),
        }
    }
}

impl<'a, R: io::Read + 'a> super::Reader<'a> for TsvReader<R> {
    type Record = TsvRecord;
    type Error = csv::Error;
    type Iterator = DeserializeRecordsIter<'a, R, TsvRecord>;

    fn records(&'a mut self) -> Self::Iterator {
        self.inner.deserialize()
    }
}

#[cfg(test)]
mod tests {
    use std::io;

    use super::super::*;

    #[test]
    fn test_tsv_records() {
        let data = b"cell	feat	dist	cell_x	cell_y	x	y
0	SCUBE3	1	475.5	630.6	13146.86026973	25793.5656964
0	SCUBE3	1	475.5	630.6	13576.7356895	38396.4273422
";
        let mut reader = tsv::TsvReader::new(io::Cursor::new(&data[..]));
        for r in reader.records() {
            match r {
                Ok(rec) => {
                    assert_eq!(rec.feature, "SCUBE3");
                    assert_eq!(rec.cell_position.x, 475.5);
                }
                Err(e) => panic!("{:?}", e),
            }
        }
    }
}