use std::io;
use std::fs;
use std::path::Path;

use csv;

use model::diffexp::CDF;


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
        writer.inner.write(["feat", measure_label, "prob"].iter()).unwrap();

        writer
    }

    pub fn write(&mut self, feature: &str, cdf: &CDF) {
        for e in cdf.iter() {
            self.inner.write([
                feature,
                &format!("{}", e.0)[..],
                &format!("{}", e.1)[..]
            ].iter()).unwrap();
        }
    }
}