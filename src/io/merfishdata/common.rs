// Copyright 2018 Johannes Köster, Till Hartmann.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bit_vec::BitVec;
use failure::Fail;
use std::iter::Iterator;
use std::path::Path;

pub type Readout = BitVec;

/// Enumeration of all accepted input formats.
pub enum Format {
    TSV,
    Binary,
    Simulation,
}

impl Format {
    /// Infer input format from given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Format {
        match path.as_ref().extension() {
            Some(e) if e == "tsv" || e == "txt" => Format::TSV,
            Some(e) if e == "sim" || e == "raw" => Format::Simulation,
            _ => Format::Binary,
        }
    }
}

pub trait MerfishRecord {
    fn cell_id(&self) -> u32;
    fn cell_name(&self) -> String;
    fn cell_pos(&self) -> (f32, f32);
    fn feature_id(&self) -> u16;
    fn feature_name(&self) -> String;
    fn hamming_dist(&self) -> u8;
    fn error_mask(&self) -> u16;
    fn codeword(&self) -> u16;
    fn readout(&self) -> u16;
    fn readout_bitvec(&self) -> Readout;
    fn count(&self) -> usize;
}

pub trait Reader<'a> {
    type Record: MerfishRecord;
    type Error: Fail;
    type Iterator: Iterator<Item=Result<Self::Record, Self::Error>> + 'a;

    fn records(&'a mut self) -> Self::Iterator;
}