use std::collections::{HashMap, HashSet};
use std::io;
use std::iter::FromIterator;
use std::path::Path;

use counter::Counter;
use failure::Fail;
use itertools::{Either, Itertools};
use rand::prelude::StdRng;
use rand::SeedableRng;
use rayon::iter::ParallelIterator;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator};

use crate::io::codebook::Codebook;
use crate::io::common::Barcode;
use crate::io::merfishdata;
use crate::io::merfishdata::tsv::{modify_record, TsvRecord};
use crate::io::merfishdata::{binary, sim, tsv, Format};
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::model::bayes::readout::Counts as MismatchCounts;
use crate::model::la::common::hamming_distance16;
use std::marker::PhantomData;
use bit_vec::BitVec;

pub struct ImageInfo {
    /// The id associated with the field-of-view in which this RNA was imaged.
    fov_id: u16,
    /// The sum of the normalized intensity associated with each pixel assigned to a given RNA.
    total_magnitude: f32,
    /// The center of the pixels associated with a given RNA in the coordinate system of each field-of-view in units of pixels.
    pixel_centroid: [u16; 2],
    /// The center of an RNA, weighted by the intensity of each pixel, in the coordinate system of each field-of-view in units of pixels.
    weighted_pixel_centroid: [f32; 2],
    /// The center of an RNA in the absolute coordinate system of the sample in microns. This quantity is calculated from the weighted_pixel_centroid.
    abs_position: [f32; 2],
    /// The number of pixels associated with an RNA.
    area: u16,
    /// The normalized intensity of each pixel for each bit in the barcode, average across all pixels assigned to an RNA.
    pixel_trace_mean: [f32; 16],
    /// The standard deviation of the normalized intensity for each bit in the barcode taken across all pixels assigned to an RNA.
    pixel_trace_std: [f32; 16],
    /// The bit at which error correction was applied. (0 if it was not applied).
    error_bit: u8,
    /// A boolean representing the type of error (0 corresponds to a '0' to '1' error; 1 corresponds to a '1' to '0' error).
    error_dir: u8,
    /// The average Euclidean distance between the pixel trace for a given RNA and the barcode to which it was matched.
    av_distance: f32,
    /// A boolean flag representing whether or not (TRUE/FALSE) an RNA was found within the nucleus of the cell.
    in_nucleus: u8,
    /// The distance (in microns) from the RNA to the closest edge of the nucleus.
    dist_nucleus: f64,
    /// The distance (in microns) from the RNA to the closest boundary of the cell.
    dist_periphery: f64,
}

pub struct CommonRecord {
    cell_name: String,
    feature_name: String,
    readout: Barcode,
    codeword: Option<Barcode>,
    count: usize,
    is_exact: bool,
    image_info: Option<ImageInfo>,
}

pub trait FromRecord<T> {
    fn from_record(other: T) -> Self;
}

// TODO impl FromRecords<BinaryRecord> for CommonRecord { ... }
impl<M: MerfishRecord> FromRecord<M> for CommonRecord {
    fn from_record(other: M) -> Self {
        CommonRecord {
            cell_name: other.cell_name(),
            feature_name: other.feature_name(),
            readout: Barcode(other.readout()),
            codeword: other.codeword().map(Barcode),
            count: other.count(),
            is_exact: other.is_exact(),
            image_info: None,
        }
    }
}

impl MerfishRecord for CommonRecord {
    fn cell_id(&self) -> u32 {
        unimplemented!()
    }

    fn cell_name(&self) -> String {
        self.cell_name.to_owned()
    }

    fn cell_pos(&self) -> (f32, f32) {
        unimplemented!()
    }

    fn feature_id(&self) -> u16 {
        unimplemented!()
    }

    fn feature_name(&self) -> String {
        self.feature_name.to_owned()
    }

    fn hamming_dist(&self) -> u8 {
        assert!(self.codeword().is_some());
        hamming_distance16(self.codeword().unwrap(), self.readout())
    }

    fn error_mask(&self) -> u16 {
        assert!(self.codeword().is_some());
        self.codeword().unwrap() ^ self.readout()
    }

    fn codeword(&self) -> Option<u16> {
        self.codeword.map(u16::from)
    }

    fn readout(&self) -> u16 {
        *self.readout
    }

    fn readout_bitvec(&self) -> Readout {
        unimplemented!()
    }

    fn count(&self) -> usize {
        self.count
    }

    fn is_exact(&self) -> bool {
        self.is_exact
    }
}

pub enum RecordIterator<R: io::Read> {
    TsvIterator { reader: tsv::Reader<R> },
    SimIterator { reader: sim::Reader<R> },
    BinaryIterator { reader: binary::Reader<R> },
}

impl RecordIterator<std::fs::File> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        let format = Format::from_path(&path);
        match format {
            Format::TSV => RecordIterator::TsvIterator {
                reader: tsv::Reader::from_file(path).unwrap(),
            },
            Format::Binary => RecordIterator::BinaryIterator {
                reader: binary::Reader::from_file(path).unwrap(),
            },
            Format::Simulation => RecordIterator::SimIterator {
                reader: sim::Reader::from_file(path).unwrap(),
            },
        }
    }
}

// TODO use
// #[enum_dispatch(MerfishRecord)]
pub enum Record {
    Tsv(tsv::TsvRecord),
    Sim(sim::SimRecord),
    Binary(binary::BinaryRecord),
}

impl MerfishRecord for Record {
    fn cell_id(&self) -> u32 {
        match self {
            Record::Tsv(r) => r.cell_id(),
            Record::Sim(r) => r.cell_id(),
            Record::Binary(r) => r.cell_id(),
        }
    }

    fn cell_name(&self) -> String {
        match self {
            Record::Tsv(r) => r.cell_name(),
            Record::Sim(r) => r.cell_name(),
            Record::Binary(r) => r.cell_name(),
        }
    }

    fn cell_pos(&self) -> (f32, f32) {
        match self {
            Record::Tsv(r) => r.cell_pos(),
            Record::Sim(r) => r.cell_pos(),
            Record::Binary(r) => r.cell_pos(),
        }
    }

    fn feature_id(&self) -> u16 {
        unimplemented!()
    }

    fn feature_name(&self) -> String {
        unimplemented!()
    }

    fn hamming_dist(&self) -> u8 {
        unimplemented!()
    }

    fn error_mask(&self) -> u16 {
        unimplemented!()
    }

    fn codeword(&self) -> Option<u16> {
        unimplemented!()
    }

    fn readout(&self) -> u16 {
        unimplemented!()
    }

    fn readout_bitvec(&self) -> BitVec<u32> {
        unimplemented!()
    }

    fn count(&self) -> usize {
        unimplemented!()
    }

    fn is_exact(&self) -> bool {
        unimplemented!()
    }
}

impl IntoIterator for RecordIterator<std::fs::File> {
    type Item = Record;
    type IntoIter = Box<dyn Iterator<Item=Record>>;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            RecordIterator::TsvIterator { mut reader } => {
                Box::new(reader.records().filter_map(Result::ok).map(Record::Tsv))
            }
            RecordIterator::SimIterator { mut reader } => {
                Box::new(reader.records().filter_map(Result::ok).map(Record::Sim))
            }
            RecordIterator::BinaryIterator { mut reader } => {
                Box::new(reader.records().filter_map(Result::ok).map(Record::Binary))
            }
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Counts {
    pub observed: Counter<(Option<Barcode>, Barcode), usize>,
    pub corrected: Counter<(Option<Barcode>, Barcode), usize>,
    pub info: Info,
}

pub(crate) type Info = HashMap<Barcode, MismatchCounts>;

pub(crate) fn into_u16(readout: &Readout) -> u16 {
    readout
        .iter()
        .rev()
        .enumerate()
        .map(|(i, bit)| (bit as u16) << i)
        .sum()
}

impl Counts {
    pub(crate) fn from_records<I: Iterator<Item=C>, C: MerfishRecord>(
        records: I,
        expressed_codewords: &[Barcode],
    ) -> HashMap<String, Self> {
        let mut cell_counts = HashMap::new();
        let mut corrections: HashMap<(Option<Barcode>, Barcode), Option<Barcode>> = HashMap::new();
        for r in records {
            let counts_in_cell = cell_counts
                .entry(r.cell_name())
                .or_insert(Counts::default());
            let (codeword, readout) = (r.codeword().map(Barcode), Barcode(r.readout()));
            *counts_in_cell
                .observed
                .entry((codeword, readout))
                .or_insert(0) += r.count();
            corrections
                .entry((codeword, readout))
                .or_insert(Self::correct_readout(
                    codeword,
                    readout,
                    &expressed_codewords,
                ));
        }
        for ((codeword, observed), corrected) in corrections {
            for (_cell, counts) in cell_counts.iter_mut() {
                let count = counts
                    .observed
                    .entry((codeword, observed))
                    .or_insert(0)
                    .clone();
                if let Some(corrected) = corrected {
                    *counts.corrected.entry((codeword, corrected)).or_insert(0) += count;
                    // The observed readout could be corrected to a codeword from the codebook
                    if observed != corrected {
                        counts
                            .info
                            .entry(corrected)
                            .or_insert(MismatchCounts {
                                exact: 0,
                                corrected: 0,
                                uncorrected: 0,
                            })
                            .corrected += count as u32;
                    } else {
                        // The observed readout is already correct
                        // (i.e. a codeword from the codebook)
                        counts
                            .info
                            .entry(observed)
                            .or_insert(MismatchCounts {
                                exact: 0,
                                corrected: 0,
                                uncorrected: 0,
                            })
                            .exact += count as u32;
                    }
                } else {
                    // The observed readout could not be corrected
                    counts
                        .info
                        .entry(observed)
                        .or_insert(MismatchCounts {
                            exact: 0,
                            corrected: 0,
                            uncorrected: 0,
                        })
                        .uncorrected += count as u32;
                }
            }
        }
        cell_counts
    }

    fn correct_readout(
        codeword: Option<Barcode>,
        readout: Barcode,
        codewords: &[Barcode],
    ) -> Option<Barcode> {
        // no correction necessary
        if codewords.contains(&readout) {
            return Some(readout);
        }

        if let Some(codeword) = codeword {
            // if codeword and readout are different (which is the case in records obtained from
            // pre-corrected sources, e.g. tsv merfish format), we already know that readout
            // should be corrected to codeword
            if codeword != readout {
                return Some(codeword);
            }
        } else {
            // if the codeword is unknown, try to find a codeword in the codebook
            // which is the only 1-neighbour of the readout

            // distances between readout and known (and expressed) codewords
            let dists = codewords
                .iter()
                .map(|&cw| (cw, hamming_distance16(*cw, *readout)))
                .collect_vec();
            // if there's exactly one expressed codeword with distance `dist`,
            // "replace" readout with that codeword.
            for dist in 1..=4 {
                let closest: Vec<Barcode> = dists
                    .iter()
                    .filter(|(_cw, d)| *d == dist)
                    .map(|(cw, _)| *cw)
                    .collect();
                if closest.len() == 1 {
                    return Some(closest[0]);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Error, Write};

    use crate::io::simple_codebook::SimpleCodebook;

    use super::*;

    #[test]
    fn test_counts() -> Result<(), Error> {
        let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        //        let path_data = "/vol/nano/merfish/la/results/obssim_A1_17_0.01_0.02.sim";
        //        let path_data = "/vol/tiny/merfish/merfishtools-evaluation/data/140genesData.1.all.txt";
        let path_codebook = "/vol/tiny/merfish/merfishtools-evaluation/codebook/140genesData.1.txt";
        let codebook = SimpleCodebook::from_file(path_codebook)?;
        let (expressed_codewords, _unexpressed_codewords) = codebook.codewords();
        let records = RecordIterator::from_path(path_data);
        let counts = Counts::from_records(records.into_iter(), &expressed_codewords);
        //        let mut file = std::fs::File::create("foo.txt")?;
        //        file.write_all(format!("{:?}", counts["0"].info).replace("}, ", "},\n").as_bytes())?;
        //        assert!(false);
        // TODO create small sim-file and codebook with known counts
        Ok(())
    }

    #[test]
    fn test_counter() {
        let a = vec![1, 2, 2, 3, 3, 3];
        let c = a.into_iter().collect::<Counter<_>>();
        let b = vec![(1, 1), (2, 1), (2, 1), (3, 3)];
        let d = b.into_iter().collect::<Counter<_>>();
        assert_eq!(&c, &d);
    }
}
