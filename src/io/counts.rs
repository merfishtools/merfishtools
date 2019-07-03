use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use std::path::Path;

use byteorder::ByteOrder;
use byteorder::NativeEndian;
use counter::Counter;
use itertools::{Either, Itertools};
use rand::prelude::StdRng;
use rand::SeedableRng;

use crate::io::codebook::Codebook;
use crate::io::merfishdata;
use crate::io::merfishdata::{binary, tsv};
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::model::la::common::hamming_distance16;

#[derive(Clone, PartialEq)]
pub(crate) struct CommonRecord {
    codeword: u16,
    readout: u16,
    cell_id: u32,
    cell_name: String,
    hamming_distance: u8,
    count: usize,
}

impl From<binary::Record> for CommonRecord {
    fn from(other: binary::Record) -> Self {
        CommonRecord {
            codeword: other.codeword(),
            readout: other.readout(),
            cell_id: other.cell_id(),
            cell_name: other.cell_name(),
            hamming_distance: other.hamming_dist(),
            count: other.count(),
        }
    }
}

impl From<tsv::Record> for CommonRecord {
    fn from(other: tsv::Record) -> Self {
        CommonRecord {
            codeword: other.codeword(),
            readout: other.readout(),
            cell_id: other.cell_id(),
            cell_name: other.cell_name(),
            hamming_distance: other.hamming_dist(),
            count: other.count(),
        }
    }
}

impl MerfishRecord for CommonRecord {
    fn cell_id(&self) -> u32 {
        self.cell_id
    }

    fn cell_name(&self) -> String {
        self.cell_name.clone()
    }

    fn cell_pos(&self) -> (f32, f32) {
        unimplemented!()
    }

    fn feature_id(&self) -> u16 {
        self.codeword
    }

    fn feature_name(&self) -> String {
        format!("{}", self.codeword)
    }

    fn hamming_dist(&self) -> u8 {
        self.hamming_distance
    }

    fn error_mask(&self) -> u16 {
        self.readout ^ self.codeword
    }

    fn codeword(&self) -> u16 {
        self.codeword
    }

    fn readout(&self) -> u16 {
        self.readout
    }

    fn readout_bitvec(&self) -> Readout {
        let mut buf = [0; 8];
        NativeEndian::write_u64(&mut buf, self.readout as u64);
        let mut readout = Readout::with_capacity(16);
        (0..16).for_each(|i| readout.push(((self.readout >> i) & 1) == 1));
        readout.truncate(16);
        if self.error_mask() > 0 {
            // error_mask == 0 <=> there is no error.
            // i.e. the actually erroneous bit is `error_bit - 1`
            let mut error_mask = Readout::with_capacity(16);
            (0..16).for_each(|i| error_mask.push(((self.error_mask() >> i) & 1) == 1));
            error_mask.truncate(16);
            // TODO swap crate bit-vec for bitvec, use xor directly.
            (0..16)
                .for_each(|i| readout.set(i, readout.get(i).unwrap() ^ error_mask.get(i).unwrap()));
        }
        readout
    }

    fn count(&self) -> usize {
        self.count
    }
}

#[derive(Debug, Clone)]
pub struct Counts {
    observed: Counter<(u16, u16), usize>,
    corrected: Counter<(u16, u16), usize>,
    info: Info,
}

pub(crate) type Info = HashMap<u16, crate::model::bayes::readout::Counts>;

pub(crate) fn into_u16(readout: &Readout) -> u16 {
    readout
        .iter()
        .rev()
        .enumerate()
        .map(|(i, bit)| (bit as u16) << i)
        .sum()
}

impl Counts {
    pub fn from_records<I: Iterator<Item=impl MerfishRecord>>(
        records: I,
        codebook: &Codebook,
    ) -> Self {
        // gather counts per readout from records
        let observed_counts = records
            .map(|r| ((r.codeword(), r.readout()), r.count()))
            .collect::<Counter<_>>();
        let (expressed_codewords, _blank_codewords): (Vec<_>, Vec<_>) = codebook
            .records()
            .iter()
            .map(|&r| (into_u16(r.codeword()), r.expressed()))
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

    pub fn from_path_grouped<P: AsRef<Path>>(
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

    pub fn from_path<P: AsRef<Path>>(path: P, codebook: &Codebook) -> Self {
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
        observed_counts: &Counter<(u16, u16), usize>,
        codewords: &HashSet<u16>,
    ) -> Vec<(u16, u16)> {
        // simple, naïve readout correction
        // If there's exactly one neighbour with hamming distance <= 4 in the codebook, use that.
        observed_counts
            .iter()
            .flat_map(|(&(codeword, readout), &count)| {
                if codeword != readout {
                    return vec![(readout, codeword); count];
                }

                // no correction necessary
                if codewords.contains(&readout) && codeword != 0 {
                    return vec![(readout, readout); count];
                }
                // if the codeword is unknown, try to find a codeword in the codebook
                // which is the only 1-neighbour of the readout
                if codeword == 0 {
                    // TODO use Option instead of 0 as a special value
                    // distances between readout and known (and expressed) codewords
                    let dists = codewords
                        .iter()
                        .map(|&cw| (cw, hamming_distance16(cw, readout)))
                        .collect_vec();
                    // if there's exactly one expressed codeword with distance `dist`,
                    // "replace" readout with that codeword.
                    for dist in 1..=2 {
                        let closest: Vec<u16> = dists
                            .iter()
                            .filter(|(_cw, d)| *d == dist)
                            .map(|(cw, _)| *cw)
                            .collect();
                        if closest.len() == 1 {
                            return vec![(readout, closest[0]); count];
                        }
                    }
                }
                vec![(readout, readout); count]
            })
            .collect_vec()
    }

    fn apply_corrections(
        observed_counts: &Counter<(u16, u16), usize>,
        corrections: &Vec<(u16, u16)>,
    ) -> (Counter<(u16, u16), usize>, Info) {
        let mut corrected_counts = observed_counts.clone();
        let mut info: Info = HashMap::new();
        for (observed_readout, corrected_readout) in corrections {
            // correct observed readout (to codeword [corrected readout]),
            // i.e. readout count is decremented
            // and codeword count incremented.
            corrected_counts
                .entry((*corrected_readout, *observed_readout))
                .and_modify(|count| *count -= 1);
            corrected_counts
                .entry((*corrected_readout, *corrected_readout))
                .and_modify(|count| *count += 1);

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

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_counts() -> Result<(), Error> {
//        let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_data = "/vol/tiny/merfish/merfishtools-evaluation/data/140genesData.1.all.txt";
        let path_codebook = "/vol/tiny/merfish/merfishtools-evaluation/codebook/140genesData.1.txt";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path(path_data, &codebook);
        dbg!(&counts);
//        assert!(false);
        Ok(())
    }

    #[test]
    fn test_grouped_counts() -> Result<(), Error> {
        // let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_data = "/vol/tiny/merfish/merfishtools-evaluation/data/140genesData.1.all.txt";
        let path_codebook = "/vol/tiny/merfish/merfishtools-evaluation/codebook/140genesData.1.txt";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path_grouped(path_data, &codebook);
        dbg!(&counts);
//        assert!(false);
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
