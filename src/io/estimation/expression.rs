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
        fs::File::create(path).map(|f| Writer::from_writer(f)).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write(["cell", "feat", "expr_ev", "expr_sd", "expr_map", "expr_ci_lower", "expr_ci_upper"].iter()).unwrap();

        writer
    }

    pub fn write(&mut self, cell: &str, feature: &str, expected_value: f64, standard_deviation: f64, map: f32, credible_interval: (f32, f32)) {
        self.inner.write([
            cell,
            feature,
            &format!("{:.*}", 2, expected_value)[..],
            &format!("{:.*}", 4, standard_deviation)[..],
            &format!("{:.*}", 0, map)[..],
            &format!("{:.*}", 2, credible_interval.0)[..],
            &format!("{:.*}", 2, credible_interval.1)[..]
        ].iter()).unwrap();
    }
}
