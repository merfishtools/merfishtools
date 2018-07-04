// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::io;
use std::fs;
use std::path::Path;
use std::collections;

use ordered_float::NotNaN;
use itertools::Itertools;
use csv;

use bio::stats::LogProb;
use bio::stats::probs::cdf;

use model::expression::{CDF, NormalizedCDF};


const HEADER: [&str; 4] = ["cell", "feat", "expr", "prob"];


/// A container for feature expression CDFs from multiple cells.
pub struct CDFs {
    inner: collections::HashMap<String, Vec<NormalizedCDF>>
}


impl CDFs {
    pub fn contains_feature(&self, feature: &str) -> bool {
        self.inner.contains_key(feature)
    }

    pub fn get(&self, feature: &str) -> Option<&Vec<NormalizedCDF>> {
        self.inner.get(feature)
    }

    pub fn get_mut(&mut self, feature: &str) -> Option<&mut Vec<NormalizedCDF>> {
        self.inner.get_mut(feature)
    }

    pub fn features(&self) -> collections::hash_map::Keys<String, Vec<NormalizedCDF>> {
        self.inner.keys()
    }

    pub fn iter(&self) -> collections::hash_map::Iter<String, Vec<NormalizedCDF>> {
        self.inner.iter()
    }
}


/// A CDF record.
#[derive(RustcDecodable, RustcEncodable)]
pub struct Record {
    pub cell: String,
    pub feature: String,
    pub expression: f64,
    pub prob: LogProb,
}


/// A reader for feature expression CDFs.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>
}


impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Option<Self> {
        fs::File::open(path).map(Reader::from_reader).expect("Error opening file.")
    }
}


impl<R: io::Read> Reader<R> {
    pub fn from_reader(rdr: R) -> Option<Self> {
        let mut inner = csv::Reader::from_reader(rdr).delimiter(b'\t');
        if inner.headers().unwrap() == HEADER {
            Some(Reader { inner })
        } else {
            None
        }
    }

    pub fn cdfs(&mut self) -> CDFs {
        let mut features = collections::HashMap::new();
        let groups = self.inner.decode().map(|res| res.ok().expect("Error reading record")).group_by(
            |rec: &Record| (rec.cell.clone(), rec.feature.clone())
        );
        for ((_, feature), records) in &groups {
            let cdf = NormalizedCDF::from_cdf(records.map(|rec| {
                cdf::Entry {
                    value: NotNaN::new(rec.expression).unwrap(),
                    prob: rec.prob,
                }
            }));
            let mut cdfs = features.entry(feature).or_insert_with(Vec::new);
            cdfs.push(cdf);
        }
        CDFs { inner: features }
    }
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        fs::File::create(path).map(Writer::from_writer).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write(HEADER.iter()).unwrap();

        writer
    }

    pub fn write(&mut self, cell: &str, feature: &str, cdf: &CDF) {
        for x in cdf.iter() {
            self.inner.write([
                cell,
                feature,
                &format!("{:.0}", x.value)[..],
                &format!("{}", *x.prob)[..]
            ].iter()).unwrap();
        }
    }
}
