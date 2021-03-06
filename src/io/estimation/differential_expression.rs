// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::fs;
use std::io;
use std::path::Path;

use bio::stats::LogProb;
use csv;

use crate::model;

/// A writer for differential expression estimates.
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
            .write_record(&[
                "feat",
                "diff_pep",
                "diff_fdr",
                "diff_bf",
                &format!("{}_map", measure_label)[..],
                &format!("{}_ci_lower", measure_label)[..],
                &format!("{}_ci_upper", measure_label)[..],
            ])
            .unwrap();

        writer
    }

    pub fn write(
        &mut self,
        feature: &str,
        differential_expression_pep: LogProb,
        fdr: LogProb,
        differential_expression_bf: model::BayesFactor,
        map: model::diffexp::DiffexpMeasure,
        credible_interval: (
            model::diffexp::DiffexpMeasure,
            model::diffexp::DiffexpMeasure,
        ),
    ) {
        self.inner
            .write_record(&[
                feature,
                &format!("{:.*e}", 2, differential_expression_pep.exp())[..],
                &format!("{:.*e}", 2, fdr.exp())[..],
                &format!("{:.*}", 2, differential_expression_bf)[..],
                &format!("{:.*}", 2, map)[..],
                &format!("{:.*}", 2, credible_interval.0)[..],
                &format!("{:.*}", 2, credible_interval.1)[..],
            ])
            .unwrap();
    }
}
