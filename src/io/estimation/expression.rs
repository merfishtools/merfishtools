// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::io;
use std::fs;
use std::path::Path;
use std::ops::Range;

use csv;


/// A writer for expression estimates.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        fs::File::create(path).map( Writer::from_writer).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::WriterBuilder::new().delimiter(b'\t').from_writer(w)
        };
        let fields = vec!["cell", "feat", "expr_map", "expr_ci_lower", "expr_ci_upper"];
        writer.inner.write_record(fields).unwrap();

        writer
    }

    pub fn write(
        &mut self,
        cell: &str,
        feature: &str,
        map: u32,
        credible_interval: &Range<&u32>
    ) {
        self.inner.write_record(&[
            cell,
            feature,
            &format!("{}", map)[..],
            &format!("{}", credible_interval.start)[..],
            &format!("{}", credible_interval.end)[..]
        ]).unwrap();
    }
}
