use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use std::path::Path;

use counter::Counter;
use itertools::{Either, Itertools};
use rand::prelude::StdRng;
use rand::SeedableRng;

use crate::io::codebook::Codebook;
use crate::io::merfishdata;
use crate::io::merfishdata::{binary, tsv};
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::io::common::Barcode;
use crate::model::la::common::hamming_distance16;

#[derive(Debug, Clone)]
pub(crate) struct Counts {
    observed: Counter<(Barcode, Barcode), usize>,
    corrected: Counter<(Barcode, Barcode), usize>,
    info: Info,
}

pub(crate) type Info = HashMap<Barcode, crate::model::bayes::readout::Counts>;

pub(crate) fn into_u16(readout: &Readout) -> u16 {
    readout
        .iter()
        .rev()
        .enumerate()
        .map(|(i, bit)| (bit as u16) << i)
        .sum()
}

impl Counts {
    pub(crate) fn from_records<I: Iterator<Item=impl MerfishRecord>>(
        records: I,
        codebook: &Codebook,
    ) -> Self {
        // gather counts per readout from records
        let observed_counts = records
            .map(|r| ((Barcode(r.codeword()), Barcode(r.readout())), r.count()))
            .collect::<Counter<_>>();
        let (expressed_codewords, _blank_codewords): (Vec<_>, Vec<_>) = codebook
            .records()
            .iter()
            .map(|&r| (Barcode(into_u16(r.codeword())), r.expressed()))
            .partition_map(|(r, expressed)| {
                if expressed {
                    Either::Left(r)
                } else {
                    Either::Right(r)
                }
            });
        let codewords = HashSet::from_iter(expressed_codewords.iter().cloned());
        let corrections = Self::corrections(&observed_counts, &codewords);
        let (corrected_counts, info) = Self::apply_corrections(&observed_counts, &corrections);
        Self {
            observed: observed_counts,
            corrected: corrected_counts,
            info,
        }
    }

    pub(crate) fn from_path_grouped<P: AsRef<Path>>(
        path: P,
        codebook: &Codebook,
    ) -> HashMap<String, Self> {
        let format = merfishdata::Format::from_path(&path);
        // TODO: remove code duplication
        match format {
            merfishdata::Format::Binary => {
                let counter = |records| Self::from_records(records, codebook);
                let grouped_records: HashMap<String, Vec<binary::Record>> =
                    merfishdata::binary::Reader::from_file(&path)
                        .unwrap()
                        .records()
                        .filter_map(Result::ok)
                        // the binary merfish format actually stores pre-corrected barcodes
                        // so we un-correct those here.
                        .map(|r| (r.cell_name(), r))
                        .into_group_map();
                grouped_records
                    .iter()
                    .map(|(cell, records)| (cell.to_owned(), counter(records.iter().cloned())))
                    .collect()
            }
            merfishdata::Format::TSV => {
                let counter = |records| Self::from_records(records, codebook);
                let mut records: Vec<tsv::Record> = merfishdata::tsv::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok)
                    .collect_vec();
                let mut rng: StdRng = StdRng::seed_from_u64(42);

                tsv::Reader::fill_in_missing_information(
                    &mut records,
                    codebook,
                    16,
                    &vec![0.04; 16],
                    &vec![0.09; 16],
                    &mut rng,
                );

                let grouped_records = records
                    .into_iter()
                    .map(|r| (r.cell_name(), r))
                    .into_group_map();
                grouped_records
                    .iter()
                    .map(|(cell, records)| (cell.to_owned(), counter(records.iter().cloned())))
                    .collect()
            }
            merfishdata::Format::Simulation => {
                let counter = |records| Self::from_records(records, codebook);
                // simulated data has not been bit corrected at all, so do nothing.
                let grouped_records = merfishdata::sim::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok)
                    .map(|r| (r.cell_name(), r))
                    .into_group_map();
                grouped_records
                    .iter()
                    .map(|(cell, records)| (cell.to_owned(), counter(records.iter().cloned())))
                    .collect()
            }
        }
    }

    pub(crate) fn from_path<P: AsRef<Path>>(path: P, codebook: &Codebook) -> Self {
        let format = merfishdata::Format::from_path(&path);
        // TODO: remove code duplication
        let counts = match format {
            merfishdata::Format::Binary => Self::from_records(
                merfishdata::binary::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                &codebook,
            ),
            merfishdata::Format::TSV => {
                let mut records: Vec<tsv::Record> = merfishdata::tsv::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok)
                    .collect_vec();
                let mut rng: StdRng = StdRng::seed_from_u64(42);

                tsv::Reader::fill_in_missing_information(
                    &mut records,
                    codebook,
                    16,
                    &vec![0.04; 16],
                    &vec![0.09; 16],
                    &mut rng,
                );
                Self::from_records(records.iter().cloned(), &codebook)
            }
            merfishdata::Format::Simulation => Self::from_records(
                merfishdata::sim::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                &codebook,
            ),
        };
        counts
    }

    fn corrections(
        observed_counts: &Counter<(Barcode, Barcode), usize>,
        codewords: &HashSet<Barcode>,
    ) -> Vec<(Barcode, Barcode)> {
        // simple, naïve readout correction
        // If there's exactly one neighbour with hamming distance <= 4 in the codebook, use that.
        observed_counts
            .iter()
            .flat_map(|(&(codeword, readout), &count)| {
                if codeword != Barcode(0) {
                    // TODO use Option instead of 0 as a special value
                    if codeword != readout {
                        return vec![(readout, codeword); count];
                    }

                    // no correction necessary
                    if codewords.contains(&readout) {
                        return vec![(readout, readout); count];
                    }
                    // if the codeword is unknown, try to find a codeword in the codebook
                    // which is the only 1-neighbour of the readout
                } else {
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
                            return vec![(readout, closest[0]); count];
                        }
                    }
                }
                vec![(readout, Barcode(0)); count]
            })
            .collect_vec()
    }

    fn apply_corrections(
        observed_counts: &Counter<(Barcode, Barcode), usize>,
        corrections: &Vec<(Barcode, Barcode)>,
    ) -> (Counter<(Barcode, Barcode), usize>, Info) {
        let mut corrected_counts = observed_counts.clone();
        let mut info: Info = HashMap::new();
        for (observed_readout, corrected_readout) in corrections {
            if *corrected_readout != Barcode(0) {
                // correct observed readout (to codeword [corrected readout]),
                // i.e. readout count is decremented
                // and codeword count incremented.
                corrected_counts
                    .entry((*corrected_readout, *observed_readout))
                    .and_modify(|count| *count -= 1);
                corrected_counts
                    .entry((*corrected_readout, *corrected_readout))
                    .and_modify(|count| *count += 1);
            }
            if observed_readout != corrected_readout {
                // if a correction has been performed
                // increment the mismatch count for the respective codeword
                info.entry(*corrected_readout)
                    .or_insert(crate::model::bayes::readout::Counts {
                        exact: 0,
                        mismatch: 0,
                    })
                    .mismatch += 1;
            } else {
                // if no correction has been performed
                // increment the exact count for the respective codeword
                // (which is identical to the readout)
                info.entry(*observed_readout)
                    .or_insert(crate::model::bayes::readout::Counts {
                        exact: 0,
                        mismatch: 0,
                    })
                    .exact += 1;
            }
        }
        (corrected_counts, info)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Error;

    use super::*;

    #[test]
    fn test_counts() -> Result<(), Error> {
        //let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_data = "/vol/nano/merfish/la/results/obssim_A1_17_0.01_0.02.sim";
//        let path_data = "/vol/tiny/merfish/merfishtools-evaluation/data/140genesData.1.all.txt";
        let path_codebook = "/vol/tiny/merfish/merfishtools-evaluation/codebook/140genesData.1.txt";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path(path_data, &codebook);
        dbg!(&counts);
        assert!(false);
        Ok(())
    }

    #[test]
    fn test_grouped_counts() -> Result<(), Error> {
        //let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_data = "/vol/tiny/merfish/merfishtools-evaluation/data/140genesData.1.all.txt";
        let path_codebook = "/vol/tiny/merfish/merfishtools-evaluation/codebook/140genesData.1.txt";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path_grouped(path_data, &codebook);
        dbg!(&counts);
        //assert!(false);
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
