use std::io;
use std::fs;
use std::path::Path;

use csv;

use bio::stats::logprobs::LogProb;

use model::foldchange::LogFC;


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
        writer.inner.write(["feat", "diff_pep", "log2fc_ev", "log2fc_sd"].iter()).unwrap();

        writer
    }

    pub fn write(&mut self, feature: &str, differential_expression_pep: LogProb, expected_value: LogFC, standard_deviation: LogFC) {
        self.inner.write([
            feature.to_owned(),
            format!("{:e}", differential_expression_pep.exp()),
            format!("{:.*}", 2, expected_value),
            format!("{:.*}", 4, standard_deviation)
        ].iter()).unwrap();
    }
}
