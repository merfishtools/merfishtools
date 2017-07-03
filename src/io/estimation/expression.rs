// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::io;
use std::fs;
use std::path::Path;
use std::ops::Range;

use ordered_float::NotNaN;
use csv;


/// A writer for expression estimates.
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    print_naive: bool
}


impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, print_naive: bool) -> Self {
        fs::File::create(path).map(|f| Writer::from_writer(f, print_naive)).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W, print_naive: bool) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t'),
            print_naive: print_naive
        };
        let mut fields = vec!["cell", "feat", "expr_map", "expr_ci_lower", "expr_ci_upper"];
        if print_naive {
            fields.push("expr_naive");
        }
        writer.inner.write(fields.iter()).unwrap();

        writer
    }

    pub fn write(
        &mut self,
        cell: &str,
        feature: &str,
        map: NotNaN<f64>,
        credible_interval: Range<&NotNaN<f64>>,
        naive_estimate: u32
    ) {
        if self.print_naive {
            self.inner.write([
                cell,
                feature,
                &format!("{}", map)[..],
                &format!("{}", credible_interval.start)[..],
                &format!("{}", credible_interval.end)[..],
                &format!("{}", naive_estimate)[..]
            ].iter()).unwrap();
        } else {
            self.inner.write([
                cell,
                feature,
                &format!("{}", map)[..],
                &format!("{}", credible_interval.start)[..],
                &format!("{}", credible_interval.end)[..]
            ].iter()).unwrap();
        }

    }
}
