use std::io;
use std::fs;
use std::path::Path;

use csv;
use itertools::Itertools;

use bio::stats::logprobs::LogProb;

use model::foldchange::{LogFC, PMF};
use model;


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
        fs::File::create(path).map(|f| Writer::from_writer(f)).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write(["feat", "log2fc", "prob"].iter()).unwrap();

        writer
    }

    pub fn write(&mut self, feature: &str, pmf: &PMF) {
        let pmf = pmf.iter().filter(|e| e.prob >= model::MIN_PROB)
                            .sorted_by(|a, b| a.value.partial_cmp(&b.value).unwrap());

        for e in pmf.iter() {
            self.inner.write([feature, &format!("{}", e.value)[..], &format!("{}", e.prob)[..]].iter()).unwrap();
        }
    }
}
