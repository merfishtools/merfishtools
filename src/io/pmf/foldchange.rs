use std::io;
use std::fs;
use std::path::Path;

use csv;

use bio::stats::logprobs::LogProb;

use model::foldchange::{LogFC, PMF};


#[derive(RustcEncodable, RustcDecodable)]
pub struct Record {
    pub feature: String,
    pub foldchange: LogFC,
    pub prob: LogProb
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        fs::File::open(path).map(|f| Writer::from_writer(f)).ok().expect("Error opening file.")
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write(["feat", "log2fc", "prob"].iter()).ok().expect("Error writing header.");

        writer
    }

    pub fn write(&mut self, feature: &str, pmf: &PMF) {
        let mut record = Record {
            feature: feature.to_owned(),
            foldchange: 0.0,
            prob: 0.0
        };

        for &(fc, prob) in pmf.iter() {
            record.foldchange = fc;
            record.prob = prob;
            self.inner.encode(&record).ok().expect("Error writing record.");
        }
    }
}
