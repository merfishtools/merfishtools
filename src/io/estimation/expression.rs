use std::io;
use std::fs;
use std::path::Path;

use csv;


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
        writer.inner.write(["expmnt", "cell", "feat", "expr_ev", "expr_sd"].iter()).ok().expect("Error writing header.");

        writer
    }

    pub fn write(&mut self, experiment: u32, cell: u32, feature: &str, expected_value: f64, standard_deviation: f64) {
        self.inner.write([
            format!("{}", experiment),
            format!("{}", cell),
            feature.to_owned(),
            format!("{:.*}", 2, expected_value),
            format!("{:.*}", 4, standard_deviation)
        ].iter()).ok().expect("Error writing record.");
    }
}
