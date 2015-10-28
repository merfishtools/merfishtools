use std::collections::HashMap;
use std::hash::Hash;
use std::collections::hash_map;
use std::io;

use itertools::Itertools;

use bio::stats::LogProb;

use io::pmf;

pub type PMF = pmf::PMF<u32>;


#[derive(RustcEncodable, RustcDecodable, PartialEq, Eq, Clone, Hash)]
pub struct Cell {
    pub experiment: u32,
    pub cell: u32
}


pub struct Features {
    inner: HashMap<String, Vec<PMF>>
}


impl Features {
    pub fn contains_feature(&self, feature: &str) -> bool {
        self.inner.contains_key(feature)
    }

    pub fn get(&self, feature: &str) -> Option<&Vec<PMF>> {
        self.inner.get(feature)
    }

    pub fn features(&self) -> hash_map::Keys<String, Vec<PMF>>  {
        self.inner.keys()
    }

    pub fn iter(&self) -> hash_map::Iter<String, Vec<PMF>> {
        self.inner.iter()
    }
}


#[derive(RustcEncodable, RustcDecodable, PartialEq, Eq, Clone, Hash)]
pub struct Record {
    pub cell: Cell,
    pub feature: String
}


pub fn features<R: io::Read>(mut reader: pmf::Reader<R>) -> Features {
    let mut features = HashMap::new();
    for (record, records) in reader.records::<Record, u32>().map(|res| res.ok().expect("Error reading record")).group_by(
        |rec| rec.feature.clone()
    ) {
        let pmf = PMF::new(records.iter().map(|rec| (rec.value.clone(), rec.prob.clone())).collect_vec());
        if !features.contains_key(&record.feature) {
            features.insert(record.feature.clone(), Vec::new());
        }
        let pmfs = features.get_mut(&record.feature).unwrap();
        pmfs.push(pmf);
    }
    Features { inner: features }
}
