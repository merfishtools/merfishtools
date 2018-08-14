// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fs;
use std::io;
use std::path::Path;

use csv;

use model::diffexp::CDF;

/// A writer for differential expression CDFs.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
}

impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, measure_label: &str) -> Self {
        fs::File::create(path)
            .map(|f| Writer::from_writer(f, measure_label))
            .unwrap()
    }
}

impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W, measure_label: &str) -> Self {
        let mut writer = Writer {
            inner: csv::WriterBuilder::new().delimiter(b'\t').from_writer(w),
        };
        writer
            .inner
            .write_record(&["feat", measure_label, "prob"])
            .unwrap();

        writer
    }

    pub fn write(&mut self, feature: &str, cdf: &CDF) {
        for e in cdf.iter() {
            self.inner
                .write_record(&[
                    feature,
                    &format!("{}", e.value)[..],
                    &format!("{}", *e.prob)[..],
                ])
                .unwrap();
        }
    }
}
