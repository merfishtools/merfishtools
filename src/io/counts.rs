use std::collections::HashMap;
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

#[derive(Debug)]
pub struct Counts {
    observed: Counter<u16, usize>,
    corrected: Counter<u16, usize>,
}

pub(crate) fn into_u16(readout: &Readout) -> u16 {
    readout
        .iter()
        .rev()
        .enumerate()
        .map(|(i, bit)| (bit as u16) << i)
        .sum()
}

impl Counts {
    pub fn from_records<I: Iterator<Item = impl MerfishRecord>>(
        records: I,
        codebook: &Codebook,
    ) -> Self {
        // gather counts per readout from records
        let observed_counts = records
            .map(|r| (r.readout(), r.count()))
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
        let corrected_counts = Self::correct_counts(&observed_counts, &expressed_codewords);
        Self {
            observed: observed_counts,
            corrected: corrected_counts,
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

                let grouped_records = records.into_iter().map(|r| (r.cell_name(), r)).into_group_map();
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
            merfishdata::Format::TSV => Self::from_records(
                merfishdata::tsv::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                &codebook,
            ),
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

    fn correct_counts(
        observed_counts: &Counter<u16, usize>,
        codewords: &[u16],
    ) -> Counter<u16, usize> {
        // simple, naïve readout correction
        // If there's exactly one neighbour with hamming distance <= 4 in the codebook, use that.
        let corrections = observed_counts
            .keys()
            .map(|&readout| {
                // no correction necessary
                if codewords.contains(&readout) {
                    return (readout, None);
                }
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
                        dbg!((&readout, &closest[0]));
                        return (readout, Some(closest[0]));
                    }
                }
                (readout, None)
            })
            .collect_vec();

        let mut corrected_counts = observed_counts.clone();
        for (readout, codeword) in corrections {
            if let Some(cw) = codeword {
                corrected_counts
                    .entry(readout)
                    .and_modify(|count| *count -= 1);
                corrected_counts.entry(cw).and_modify(|count| *count += 1);
            };
        }
        corrected_counts
    }
}

#[cfg(test)]
mod tests {
    use std::io::Error;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_counts() -> Result<(), Error> {
        let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_codebook = "/home/hartmann/projects/merfresh/codebook.tsv";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path(path_data, &codebook);
        dbg!(&counts);
        assert!(false);
        Ok(())
    }

    #[test]
    fn test_grouped_counts() -> Result<(), Error> {
        let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_codebook = "/home/hartmann/projects/merfresh/codebook.tsv";
        let codebook = Codebook::from_file(path_codebook)?;
        let counts = Counts::from_path_grouped(path_data, &codebook);
        dbg!(&counts);
        assert!(false);
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
