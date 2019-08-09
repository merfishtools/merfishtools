use std::collections::{HashMap, HashSet};
use std::io;
use std::iter::FromIterator;
use std::marker::PhantomData;
use std::path::Path;

use bit_vec::BitVec;
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
use crate::io::merfishdata::common::{Record, RecordReader};
use crate::io::merfishdata::tsv::{modify_record, TsvRecord};
use crate::io::merfishdata::{binary, sim, tsv, Format};
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::model::bayes::readout::Counts as MismatchCounts;
use crate::model::la::common::hamming_distance16;

impl RecordReader<std::fs::File> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Self {
        let format = Format::from_path(&path);
        match format {
            Format::TSV => RecordReader::TsvIterator {
                reader: tsv::TsvReader::from_file(path).unwrap(),
            },
            Format::Binary => RecordReader::BinaryIterator {
                reader: binary::BinaryReader::from_file(path).unwrap(),
            },
            Format::Simulation => RecordReader::SimIterator {
                reader: sim::SimReader::from_file(path).unwrap(),
            },
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
    pub(crate) fn from_records<I: Iterator<Item = C>, C: MerfishRecord>(
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
        let records = RecordReader::from_path(path_data);
        let counts = Counts::from_records(records.into_iter(), &expressed_codewords);
        let mut file = std::fs::File::create("foo.txt")?;
        file.write_all(
            format!("{:?}", counts["0"].info)
                .replace("}, ", "},\n")
                .as_bytes(),
        )?;
        assert!(false);
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
