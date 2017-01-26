// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::io;
use std::fs;
use std::path::Path;
use std::collections::{HashMap, hash_map};

use csv;
use itertools::Itertools;
use bio::alignment::distance;

/// A codebook record.
#[derive(RustcDecodable)]
pub struct Record {
    pub feature: String,
    pub codeword: String,
    pub expressed: u8
}


/// Codebook representation.
pub struct Codebook {
    inner: HashMap<String, [u32; 4]>
}


impl Codebook {
    pub fn neighbors(&self, feature: &str, dist: u8) -> u32 {
        self.inner.get(feature)
                  .expect(&format!("Error: Feature {} not in codebook.", feature))[dist as usize - 1]
    }

    pub fn features(&self) -> hash_map::Keys<String, [u32; 4]> {
        self.inner.keys()
    }

    pub fn contains(&self, feature: &str) -> bool {
        self.inner.contains_key(feature)
    }
}


/// Reader for codebook definitions.
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    dist: u8
}


impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, dist: u8) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::from_reader(f, dist))
    }
}


impl<R: io::Read> Reader<R> {
    pub fn from_reader(rdr: R, dist: u8) -> Self {
        Reader {
            inner: csv::Reader::from_reader(rdr).delimiter(b'\t'),
            dist: dist
        }
    }

    pub fn codebook(&mut self) -> Codebook {
        let records = self.inner.decode::<Record>().map(|record| record.unwrap()).collect_vec();
        let mut inner = HashMap::new();
        for record in records.iter() {
            if record.expressed == 1 {
                inner.insert(record.feature.clone(), [0; 4]);
            }
        }
        for a in records.iter() {
            let codeword = a.codeword.as_bytes();
            let neighbors = inner.get_mut(&a.feature).unwrap();
            for b in records.iter() {
                if b.feature != a.feature {
                    let dist = distance::hamming(codeword, b.codeword.as_bytes()).unwrap();
                    assert!(dist >= self.dist as u32, "Unexpected hamming distance {} (>={} allowed).", dist, self.dist);
                    if dist <= 4 && b.expressed == 1 {
                        neighbors[(dist - 1) as usize] += 1;
                    }
                }
            }
        }
        Codebook { inner: inner }
    }
}
