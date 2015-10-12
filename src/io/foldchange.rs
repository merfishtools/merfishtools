use std::io;
use std::fs;
use std::path::Path;

use csv;
use itertools;

use bio::stats::logprobs::{LogProb, Prob};


#[derive(RustcEncodable)]
pub struct Record {
    pub rna: String,
    pub map: f64,
    pub pep: Prob,
    pub cdf: Vec<Prob>
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    min_fs: f64,
    max_fc: f64
}


impl Writer<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, max_x: u32) -> io::Result<Self> {
        fs::File::open(path).map(|f| Writer::from_writer(f, max_x))
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W, min_fc: f64, max_fc: f64) -> Self {
        Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t'),
            min_fc: min_fc,
            max_fc: max_fc
        }
    }

    pub fn write_header(&mut self) -> csv::Result<()> {
        let mut hdr = vec!["RNA".to_string()];
        hdr.extend(itertools::linspace(self.min_fc, self.max_fc, 10).map(|fc| format!("{}", fc)));
        self.inner.write(hdr.into_iter())
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
