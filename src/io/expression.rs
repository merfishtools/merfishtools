use std::io;
use std::fs;
use std::path::Path;

use csv;

use bio::stats::logprobs::Prob;


#[derive(RustcEncodable)]
pub struct Record {
    pub experiment: u32,
    pub cell: u32,
    pub feature: String,
    pub expr: u32,
    pub prob: Prob
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
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
        self.inner.write(["Expmt", "Cell", "Feat", "Expr", "Prob"].into_iter())
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
