use std::io;
use std::fs;
use std::path::Path;

use csv;

use bio::stats::logprobs::LogProb;


#[derive(RustcEncodable)]
pub struct Record {
    pub experiment: u32,
    pub cell: u32,
    pub rna: String,
    pub pmf: Vec<Prob>
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    max_x: u32
}


impl Writer<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, max_x: u32) -> io::Result<Self> {
        fs::File::open(path).map(|f| Writer::from_writer(f, max_x))
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W, max_x: u32) -> Self {
        Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t'),
            max_x: max_x
        }
    }

    pub fn write_header(&mut self) -> csv::Result<()> {
        let mut hdr = vec!["Experiment".to_string(), "Cell".to_string(), "RNA".to_string()];
        hdr.extend((0..self.max_x).map(|x| format!("{}", x)));
        self.inner.write(hdr.into_iter())
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
