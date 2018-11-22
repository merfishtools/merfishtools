// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use num::CheckedAdd;
use std::collections::{hash_map, HashMap, HashSet};
use std::path::Path;

use bit_vec::BitVec;
use csv;
use itertools::Itertools;
use petgraph::prelude::*;

pub type Codeword = BitVec<u32>;

pub fn parse_codeword(codeword: &[u8]) -> Codeword {
    let mut _codeword = BitVec::with_capacity(codeword.len());
    for &b in codeword {
        if b == b'1' {
            _codeword.push(true);
        } else if b == b'0' {
            _codeword.push(false)
        } else {
            panic!("invalid codeword {:?}", codeword)
        }
    }

    _codeword
}

/// A codebook record.
#[derive(Clone, Debug)]
pub struct Record {
    name: String,
    codeword: BitVec<u32>,
    expressed: bool,
}

impl Record {
    pub fn noise(len: usize) -> Self {
        Record {
            name: "noise".to_owned(),
            codeword: BitVec::from_elem(len, false),
            expressed: true, // there is always some noise
        }
    }

    /// Create new record.
    pub fn new(name: String, codeword: &[u8], expressed: bool) -> Self {
        Record {
            name,
            codeword: parse_codeword(codeword),
            expressed,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn codeword(&self) -> &BitVec<u32> {
        &self.codeword
    }

    pub fn expressed(&self) -> bool {
        self.expressed
    }

    /// Get distance to other codeword.
    pub fn dist(&self, other: &Record) -> u8 {
        let mut dist = 0;
        for (a, b) in self.codeword.blocks().zip(other.codeword.blocks()) {
            dist = dist
                .checked_add(&((a ^ b).count_ones() as u8))
                .expect("bug: overflow when calculating hamming distance");
        }
        dist
    }

    pub fn diff(&self, other: &Record) -> BitVec {
        let mut diff = BitVec::new();
        {
            let s = unsafe { diff.storage_mut() };
            for (a, b) in self.codeword.blocks().zip(other.codeword.blocks()) {
                s.push(a ^ b);
            }
        }
        unsafe {
            diff.set_len(self.codeword.len());
        }

        diff
    }
}

pub type FeatureID = usize;

/// Codebook representation.\
#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct Codebook {
    index: HashMap<String, FeatureID>,
    idmap: HashMap<u16, String>,
    graph: UnGraph<Record, ()>,
    pub min_dist: u8,
    // Number of bits in the codewords.
    pub N: u8,
    // Number of 1-bits in the codewords.
    pub m: u8,
    noise_record: Record,
}

impl Codebook {
    /// Read from a given file path.
    #[allow(non_snake_case)]
    pub fn from_file<P: AsRef<Path>>(path: P) -> csv::Result<Self> {
        let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;

        let mut N = None;
        let mut m = None;
        let mut graph = Graph::default();
        let mut idmap = HashMap::new();
        let index = {
            let mut index = HashMap::new();
            for (i, rec) in rdr.deserialize().enumerate() {
                let (feature, codeword, expressed): (String, String, u8) = rec?;
                idmap.insert((i + 1) as u16, feature.clone());
                let expressed = expressed == 1;

                let rec = Record::new(feature.clone(), codeword.as_bytes(), expressed);
                let n = rec.codeword.len();
                let m_ = rec.codeword.iter().filter(|b| *b).count();

                if let Some(N) = N {
                    assert!(n == N, "unexpected codeword length: {} vs {}", n, N);
                } else {
                    N = Some(n);
                }
                if let Some(m) = m {
                    assert!(m == m_, "unexpected number of 1-bits: {} vs {}", m_, m);
                } else {
                    m = Some(m_);
                }

                let idx = graph.add_node(rec);
                index.insert(feature, idx.index());
            }
            index
        };

        let dist = |a, b, graph: &UnGraph<Record, ()>| {
            let rec_a = graph.node_weight(NodeIndex::new(a)).unwrap();
            let rec_b = graph.node_weight(NodeIndex::new(b)).unwrap();
            rec_a.dist(&rec_b)
        };

        let min_dist = index
            .values()
            .tuple_combinations()
            .map(|(&a, &b)| dist(a, b, &graph))
            .min()
            .unwrap();

        for (&a, &b) in index.values().tuple_combinations() {
            let d = dist(a, b, &graph);
            if d == min_dist {
                graph.add_edge(NodeIndex::new(a), NodeIndex::new(b), ());
            }
        }
        Ok(Codebook {
            index,
            idmap,
            graph,
            min_dist,
            N: N.unwrap() as u8,
            m: m.unwrap() as u8,
            noise_record: Record::noise(N.unwrap()),
        })
    }

    pub fn noise(&self) -> &Record {
        &self.noise_record
    }

    pub fn record(&self, feature: FeatureID) -> &Record {
        self.graph.node_weight(NodeIndex::new(feature)).unwrap()
    }

    pub fn get_name(&self, feature_id: u16) -> Option<String> {
        self.idmap.get(&feature_id).cloned()
    }

    pub fn get_id(&self, feature: &str) -> FeatureID {
        *self
            .index
            .get(feature)
            .expect("bug: feature not in codebook")
    }

    pub fn features(&self) -> hash_map::Keys<String, FeatureID> {
        self.index.keys()
    }

    pub fn records<'a>(&'a self) -> Vec<&'a Record> {
        //self.graph.node_weights_mut()
        self.graph
            .node_indices()
            .map(|n| self.graph.node_weight(n).unwrap())
            .collect_vec()
    }

    pub fn contains(&self, feature: &str) -> bool {
        self.index.contains_key(feature)
    }

    /// Return neighbors with given shortest distance to feature.
    /// If a neighbor is e.g. distance 4 and distance 2, and we require distance 4, it will not be
    /// returned because the shortest distance is 2.
    pub fn neighbors(&self, feature: FeatureID, dist: u8) -> Vec<FeatureID> {
        assert!(
            dist % self.min_dist == 0,
            "unsupported distance: only multiples of {} allowed",
            self.min_dist
        );
        assert!(dist == 4 || dist == 2, "unsupported distance {}", dist);
        let feature = NodeIndex::new(feature);

        if dist == self.min_dist {
            self.graph
                .neighbors(feature)
                .map(|n| n.index())
                .collect_vec()
        } else if dist == self.min_dist * 2 {
            let mut neighbors = HashSet::new();
            for n in self.graph.neighbors(feature) {
                for n_ in self.graph.neighbors(n) {
                    if self.graph.find_edge(feature, n_).is_none() && n_ != feature {
                        // consider if not a direct neighbor nor feature itself
                        neighbors.insert(n_.index());
                    }
                }
            }
            neighbors.into_iter().collect_vec()
        } else {
            panic!(
                "unsupported distance {}, only first and second order neighbors supported",
                dist
            );
        }
    }

    pub fn len(&self) -> usize {
        self.index.len()
    }

    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }
}
