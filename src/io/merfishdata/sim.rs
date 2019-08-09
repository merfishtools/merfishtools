use std::fs;
use std::io;
use std::iter::Map;
use std::path::Path;

use crate::io::counts::{CommonRecord, FromRecord};
use crate::io::merfishdata::{MerfishRecord, Readout};
use crate::model::la::common::hamming_distance16;
use crate::simulation::binary;
use csv::DeserializeRecordsIter;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SimRecord {
    pub cell: usize,
    #[serde(with = "binary", rename = "barcode")]
    pub readout: u16,
    pub count: usize,
    #[serde(default, skip)]
    pub codeword: Option<u16>,
}

impl MerfishRecord for SimRecord {
    fn cell_id(&self) -> u32 {
        self.cell as u32
    }

    fn cell_name(&self) -> String {
        format!("{}", self.cell)
    }

    fn cell_pos(&self) -> (f32, f32) {
        unimplemented!()
    }

    fn feature_id(&self) -> u16 {
        self.readout
    }

    fn feature_name(&self) -> String {
        format!("{}", self.readout)
    }

    fn hamming_dist(&self) -> u8 {
        // FIXME
        if let Some(codeword) = self.codeword {
            hamming_distance16(codeword, self.readout)
        } else {
            let d = self.readout.count_ones() as isize - 4;
            if d < 0 {
                (-d) as u8
            } else {
                d as u8
            }
        }
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
        self.count
    }

    fn is_exact(&self) -> bool {
        if self.codeword.is_some() {
            self.error_mask() == 0
        } else {
            false
        }
    }
}

/// A reader for MERFISH raw data.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
}

impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    pub fn new(rdr: R) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new().delimiter(b'\t').from_reader(rdr),
        }
    }
}

impl<'a, R: io::Read + 'a> super::Reader<'a> for Reader<R> {
    type Record = SimRecord;
    type Error = csv::Error;
    type Iterator = DeserializeRecordsIter<'a, R, SimRecord>;

    fn records(&'a mut self) -> Self::Iterator {
        self.inner.deserialize()
    }
}
