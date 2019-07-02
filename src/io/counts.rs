use std::path::Path;

use byteorder::ByteOrder;
use byteorder::NativeEndian;
use counter::Counter;
use itertools::{Either, Itertools};
use rand::distributions::WeightedIndex;
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

use crate::io::codebook::Codebook;
use crate::io::merfishdata;
use crate::io::merfishdata::{binary, tsv};
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::model::la::common::hamming_distance16;

struct CommonRecord {
    codeword: u16,
    readout: u16,
    cell_id: u32,
    cell_name: String,
    hamming_distance: u8,
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

    fn barcode(&self, _codebook: Option<&Codebook>) -> u16 {
        self.readout
    }

    fn readout(&self) -> Readout {
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
        1
    }
}

#[derive(Debug)]
pub struct Counts {
    observed: Counter<u16, usize>,
    corrected: Counter<u16, usize>,
}

fn into_u16(readout: &Readout) -> u16 {
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
            .map(|r| (r.barcode(Some(&codebook)), r.count()))
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

    pub fn from_path<P: AsRef<Path>>(path: P, codebook: Codebook) -> Self {
        let format = merfishdata::Format::from_path(&path);
        // TODO: remove code duplication
        let counts = match format {
            merfishdata::Format::Binary => Self::from_records(
                merfishdata::binary::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok)
                    // the binary merfish format actually stores pre-corrected barcodes
                    // so we un-correct those here.
                    .map(Self::uncorrect_binary_format_record),
                &codebook,
            ),
            merfishdata::Format::TSV => {
                let mut rng = StdRng::seed_from_u64(42);
                let uncorrector = |r: tsv::Record| {
                    Self::uncorrect_tsv_format_record(
                        r,
                        &codebook,
                        16,
                        &vec![0.04; 16],
                        &vec![0.09; 16],
                        &mut rng,
                    )
                };
                Self::from_records(
                    // The old csv merfish format only has corrected readouts.
                    // Since it does not say which bit has been corrected,
                    // un-correcting requires random correction according to
                    // positional error probabilities.
                    merfishdata::tsv::Reader::from_file(&path)
                        .unwrap()
                        .records()
                        .filter_map(Result::ok)
                        .map(uncorrector),
                    &codebook,
                )
            }
            merfishdata::Format::Simulation => Self::from_records(
                // simulated data has not been bit corrected at all, so do nothing.
                merfishdata::sim::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                &codebook,
            ),
        };
        counts
    }

    fn uncorrect_binary_format_record(record: binary::Record) -> binary::Record {
        let mut uncorrected = record.clone();
        uncorrected.barcode = if record.hamming_dist() == 1 {
            record.barcode ^ record.error_mask() as u64
        } else {
            record.barcode
        };
        uncorrected
    }

    fn uncorrect_tsv_format_record(
        record: tsv::Record,
        codebook: &Codebook,
        num_bits: usize,
        p0: &[f32],
        p1: &[f32],
        rng: &mut StdRng,
    ) -> CommonRecord {
        if record.hamming_dist() == 1 {
            let codeword = record.barcode(Some(codebook));
            let weights: Vec<_> = (0..num_bits)
                .map(|i| {
                    let bit = (codeword >> i) & 1;
                    match bit {
                        0 => p0[i],
                        1 => p1[i],
                        _ => panic!("bug: bit can only be 0 or 1"),
                    }
                })
                .collect();
            let wc = WeightedIndex::new(weights).unwrap();
            let error_bit = rng.sample(&wc) as u8;
            CommonRecord {
                codeword,
                readout: codeword ^ (1 << error_bit as u16),
                cell_id: record.cell_id(),
                cell_name: record.cell_name(),
                hamming_distance: 1,
            }
        } else if record.hamming_dist() == 0 {
            CommonRecord {
                codeword: record.barcode(Some(codebook)),
                readout: record.barcode(Some(codebook)),
                cell_id: record.cell_id(),
                cell_name: record.cell_name(),
                hamming_distance: 0,
            }
        } else {
            panic!(
                "Unsupported number of erroneous bits {}",
                record.hamming_dist()
            );
        }
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
                for dist in 1..=4 {
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
    fn test_add() -> Result<(), Error> {
        let path_data = "/vol/tiny/merfish/raw/rep2/assigned_blist.bin";
        let path_codebook = "/home/hartmann/projects/merfresh/codebook.tsv";
        let counts = Counts::from_path(path_data, Codebook::from_file(path_codebook)?);
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
