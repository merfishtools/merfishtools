use std::io;
use std::fs;
use std::path::Path;

use csv;
use itertools;

use bio::stats::logprobs::{LogProb, Prob};


#[derive(RustcEncodable)]
pub struct Record {
    pub feature: String,
    pub fc: f64,
    pub prob: Prob
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, verbose: bool) -> io::Result<Self> {
        fs::File::open(path).map(|f| Writer::from_writer(f))
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        }
    }

    pub fn write_header(&mut self) -> csv::Result<()> {
        self.inner.write(["Feat", "logFC", "Prob"].into_iter())
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
