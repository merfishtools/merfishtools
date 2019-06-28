use crate::io::codebook::Codebook;
use crate::io::merfishdata;
use crate::io::merfishdata::{MerfishRecord, Reader};
use crate::model::la::common::hamming_distance;
use itertools::Itertools;
use std::collections::HashMap;
use std::io;
use std::path::Path;

pub struct Counts {
    observed: HashMap<u16, usize>,
    corrected: HashMap<u16, usize>,
}

impl Counts {
    pub fn from_path<P: AsRef<Path>>(path: P, codebook: Codebook) -> Self {
        let format = merfishdata::Format::from_path(&path);
        // TODO: remove code duplication
        let mut records: HashMap<u16, Vec<u8>> = match format {
            merfishdata::Format::Binary => merfishdata::binary::Reader::from_file(&path)
                .unwrap()
                .records()
                .into_iter()
                .filter_map(Result::ok)
                .map(|r| (r.barcode(Some(&codebook)), r.hamming_dist()))
                .into_group_map(),
            merfishdata::Format::TSV => merfishdata::tsv::Reader::from_file(&path)
                .unwrap()
                .records()
                .into_iter()
                .filter_map(Result::ok)
                .map(|r| (r.barcode(Some(&codebook)), r.hamming_dist()))
                .into_group_map(),
            merfishdata::Format::Simulation => merfishdata::sim::Reader::from_file(&path)
                .unwrap()
                .records()
                .into_iter()
                .filter_map(Result::ok)
                .map(|r| (r.barcode(Some(&codebook)), r.hamming_dist()))
                .into_group_map(),
        };
        let observed: HashMap<u16, usize> = records
            .iter()
            .map(|(readout, dist)| (*readout, dist.len()))
            .collect();
        let codewords: Vec<(u16, bool)> = codebook.records().iter().map(|&r| {
            (r.codeword()
                 .iter()
                 .rev()
                 .enumerate()
                 .map(|(i, bit)| (bit as u16) << i)
                 .sum(),
             r.expressed())
        }).collect();
        let expressed = codewords.iter().filter(|(_, b)| *b).map(|(r, _)| *r).collect_vec();
        for cw in expressed {
            for (barcode, _) in observed {  
                let dist = hamming_distance(cw as usize, barcode as usize);
            }
        }

        unimplemented!()
    }
}
