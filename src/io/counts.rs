use crate::io::codebook::Codebook;
use crate::io::merfishdata;
use crate::io::merfishdata::binary;
use crate::io::merfishdata::{MerfishRecord, Reader, Readout};
use crate::model::la::common::hamming_distance;
use counter::Counter;
use itertools::{Either, Itertools};
use std::path::Path;

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
        codebook: Codebook,
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
                    .map(|r: binary::Record| {
                        let mut uncorrected = r.clone();
                        uncorrected.barcode = if r.hamming_dist() == 1 {
                            r.barcode ^ (1 << r.error_bit)
                        } else {
                            r.barcode
                        };
                        uncorrected
                    }),
                codebook,
            ),
            merfishdata::Format::TSV => Self::from_records(
                merfishdata::tsv::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                codebook,
            ),
            merfishdata::Format::Simulation => Self::from_records(
                merfishdata::sim::Reader::from_file(&path)
                    .unwrap()
                    .records()
                    .filter_map(Result::ok),
                codebook,
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
                    .map(|&cw| (cw, hamming_distance(cw as usize, readout as usize)))
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
