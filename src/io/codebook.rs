// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::path::Path;
use std::collections::{HashMap, hash_map, HashSet};
use num::CheckedAdd;

use csv;
use itertools::Itertools;
use bit_vec::BitVec;
use petgraph::prelude::*;


pub type Codeword = BitVec<u32>;


/// A codebook record.
#[derive(Debug)]
pub struct Record {
    name: String,
    codeword: BitVec<u32>
}


impl Record {
    /// Create new record.
    pub fn new(name: String, codeword: &[u8]) -> Self {
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

        Record {
            name: name,
            codeword: _codeword
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn codeword(&self) -> &BitVec<u32> {
        &self.codeword
    }

    /// Get distance to other codeword.
    pub fn dist(&self, other: &Record) -> u8 {
        let mut dist = 0;
        for (a, b) in self.codeword.blocks().zip(other.codeword.blocks()) {
            dist = dist.checked_add(&((a ^ b).count_ones() as u8)).expect(
                "bug: overflow when calculating hamming distance"
            );
        }
        dist
    }
}


/// Codebook representation.\
#[allow(non_snake_case)]
pub struct Codebook {
    index: HashMap<String, NodeIndex<u32>>,
    graph: UnGraph<Record, ()>,
    pub min_dist: u8,
    // Number of bits in the codewords.
    pub N: u8,
    // Number of 1-bits in the codewords.
    pub m: u8
}



impl Codebook {
    /// Read from a given file path.
    #[allow(non_snake_case)]
    pub fn from_file<P: AsRef<Path>>(path: P) -> csv::Result<Self> {
        let mut rdr = (csv::Reader::from_file(path)?).delimiter(b'\t');

        let mut N = None;
        let mut m = None;
        let mut graph = Graph::default();
        let index = {
            let mut index = HashMap::new();
            for rec in rdr.decode() {
                let (feature, codeword): (String, String) = rec?;

                let rec = Record::new(feature.clone(), codeword.as_bytes());
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
                index.insert(
                    feature,
                    idx
                );
            }
            index
        };

        let dist = |a, b, graph: &UnGraph<Record, ()>| {
            let rec_a = graph.node_weight(a).unwrap();
            let rec_b = graph.node_weight(b).unwrap();
            rec_a.dist(&rec_b)
        };

        let min_dist = index.values().tuple_combinations().map(|(&a, &b)| dist(a, b, &graph)).min().unwrap();

        for (&a, &b) in index.values().tuple_combinations() {
            let d = dist(a, b, &graph);
            if d == min_dist {
                graph.add_edge(a, b, ());
            }
        }

        Ok(Codebook {
            graph: graph,
            index: index,
            min_dist: min_dist,
            N: N.unwrap() as u8,
            m: m.unwrap() as u8
        })
    }

    pub fn get(&self, feature: &str) -> &Record {
        self.graph.node_weight(
            *self.index.get(feature).expect("bug: feature not in codebook")
        ).unwrap()
    }

    pub fn features(&self) -> hash_map::Keys<String, NodeIndex<u32>> {
        self.index.keys()
    }

    pub fn records(&self) -> impl Iterator<&Record> {
        self.graph.node_references().map(|n| n.weight())
    }

    pub fn contains(&self, feature: &str) -> bool {
        self.index.contains_key(feature)
    }

    pub fn neighbors(&self, feature: &str, dist: u8) -> Vec<&Record> {
        assert!(
            dist % self.min_dist == 0,
            "unsupported distance: only multiples of {} allowed",
            self.min_dist
        );

        let feature = *self.index.get(feature).unwrap();

        let mut visited = HashSet::new();
        let mut neighbors = Vec::new();
        self.dfs(feature, 0, &mut visited, &mut neighbors, dist);
        println!("{:?}", neighbors);

        neighbors
    }

    fn dfs<'a>(
        &'a self,
        node: NodeIndex<u32>,
        d: u8,
        visited: &mut HashSet<NodeIndex<u32>>,
        neighbors: &mut Vec<&'a Record>,
        dist: u8) {
        if visited.contains(&node) {
            return;
        }
        visited.insert(node);

        if d == dist {
            neighbors.push(self.graph.node_weight(node).unwrap());
            return;
        }

        for n in self.graph.neighbors(node) {
            self.dfs(n, d + self.min_dist, visited, neighbors, dist);
        }
    }
}
