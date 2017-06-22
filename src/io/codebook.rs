// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::io;
use std::fs;
use std::path::Path;
use std::collections::{HashMap, hash_map, HashSet};
use num::CheckedAdd;

use csv;
use itertools::Itertools;
use bit_vec::BitVec;
use petgraph::prelude::*;
use petgraph::visit;
use bio::alignment::distance;


#[derive(Debug)]
pub struct Record {
    pub name: String,
    pub codeword: BitVec<u32>
}


impl Record {
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


/// Codebook representation.
pub struct Codebook {
    index: HashMap<String, NodeIndex<u32>>,
    graph: UnGraph<Record, ()>,
    min_dist: u8
}



impl Codebook {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P, min_dist: u8) -> csv::Result<Self> {
        let mut rdr = (csv::Reader::from_file(path)?).delimiter(b'\t');

        let mut graph = Graph::default();
        let index = {
            let mut index = HashMap::new();
            for rec in rdr.decode() {
                let (feature, codeword): (String, String) = rec?;
                let idx = graph.add_node(Record::new(feature.clone(), codeword.as_bytes()));
                index.insert(
                    feature,
                    idx
                );
            }
            index
        };

        for &a in index.values() {
            for &b in index.values() {
                if a != b {
                    let d = {
                        let rec_a = graph.node_weight(a).unwrap();
                        let rec_b = graph.node_weight(b).unwrap();
                        rec_a.dist(&rec_b)
                    };

                    if d == min_dist {
                        graph.add_edge(a, b, ());
                    } else if d < min_dist {
                        panic!("unexpected hamming distance {} (given minimum {})", d, min_dist);
                    }
                }
            }
        }

        Ok(Codebook {
            graph: graph,
            index: index,
            min_dist: min_dist
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
