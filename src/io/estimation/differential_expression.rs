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
    pub fn from_file<P: AsRef<Path>>(path: P, measure_label: &str) -> Self {
        fs::File::create(path).map(|f| Writer::from_writer(f, measure_label)).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W, measure_label: &str) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write([
            "feat",
            "diff_pep",
            "diff_fdr",
            "diff_bf",
            &format!("{}_ev", measure_label)[..],
            &format!("{}_sd", measure_label)[..],
            &format!("{}_ci_lower", measure_label)[..],
            &format!("{}_ci_upper", measure_label)[..]
        ].iter()).unwrap();

        writer
    }

    pub fn write(
        &mut self,
        feature: &str,
        differential_expression_pep: LogProb,
        fdr: LogProb,
        differential_expression_bf: f64,
        expected_value: LogFC,
        standard_deviation: LogFC,
        credible_interval: (LogFC, LogFC)
    ) {
        self.inner.write([
            feature,
            &format!("{:.*e}", 2, differential_expression_pep.exp())[..],
            &format!("{:.*e}", 2, fdr.exp())[..],
            &format!("{:.*}", 2, differential_expression_bf)[..],
            &format!("{:.*}", 2, expected_value)[..],
            &format!("{:.*}", 4, standard_deviation)[..],
            &format!("{:.*}", 2, credible_interval.0)[..],
            &format!("{:.*}", 2, credible_interval.1)[..]
        ].iter()).unwrap();
    }
}
