use std::fs;
use std::io;
use std::path::Path;

use crate::io::merfishdata::{MerfishRecord, Readout};
use crate::simulation::binary;

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Record {
    pub cell: usize,
    #[serde(with = "binary", rename = "barcode")]
    pub readout: u16,
    pub count: usize,
    #[serde(with = "binary", skip)]
    pub codeword: u16,
}

impl MerfishRecord for Record {
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
        let d = self.readout.count_ones() as isize - 4;
        if d < 0 {
            (-d) as u8
        } else {
            d as u8
        }
    }

    fn error_mask(&self) -> u16 {
        0b0000_0000_0000_0000
    }

    fn codeword(&self) -> u16 {
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
    type Record = Record;
    type Error = csv::Error;
    type Iterator = csv::DeserializeRecordsIter<'a, R, Record>;

    fn records(&'a mut self) -> csv::DeserializeRecordsIter<'a, R, Record> {
        self.inner.deserialize()
    }
}