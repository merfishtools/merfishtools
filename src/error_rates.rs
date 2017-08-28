use std::mem;
use std::f64;

use itertools::Itertools;
use ndarray::prelude::*;
use bit_vec::BitVec;
use ordered_float::NotNaN;
use ndarray_rand::RandomExt;
use rand::distributions::Range;

use bio::stats::{Prob, LogProb};

use io::codebook::Codebook;

pub type Readout = BitVec;


/// Estimate position-wise error rates as presented by Chen et al. 2015.
pub fn estimate<I: Iterator<Item=(String, Readout)>>(
    codebook: &Codebook,
    readouts: I
) -> (Vec<Prob>, Vec<Prob>) {
    let max_diff = 0.001;
    let shape = (codebook.len(), codebook.N as usize);

    let mut observed_01_errors = Array2::from_elem(shape, 0.0);
    let mut observed_10_errors = Array2::from_elem(shape, 0.0);
    let mut exact_count = Array2::from_elem((codebook.len(), 1), 0.0);

    for (feat, readout) in readouts {
        let feature_id = codebook.get_id(&feat);
        let rec = codebook.record(feature_id);
        if !rec.expressed() {
            continue;
        }
        let encoding = rec.codeword();

        let mut err_count = 0;
        for k in 0..codebook.N as usize {
            match (encoding[k], readout[k]) {
                (false, true) => {
                    observed_01_errors[(feature_id, k)] += 1.0;
                    err_count += 1;
                },
                (true, false) => {
                    observed_10_errors[(feature_id, k)] += 1.0;
                    err_count += 1;
                },
                _ => ()
            };
        }
        assert!(
            err_count <= 1,
            "unexpected number of errors in readout: {} (<=1 expected)",
            err_count
        );
        if err_count == 0 {
            exact_count[(feature_id, 0)] += 1.0;
        }
    }

    // calculate observed (joint) error rates
    let observed_p0 = observed_01_errors / &exact_count;
    let observed_p1 = observed_10_errors / &exact_count;

    // estimate error rates
    let feat_p0 = observed_p0.map(|&o| o / (1.0 + o));
    let feat_p1 = observed_p1.map(|&o| o / (1.0 + o));

    let mut p0 = Vec::new();
    let mut p1 = Vec::new();
    for k in 0..codebook.N as usize {
        let p = |feat_p: &Array2<f64>, bit: bool| {
            let mut count = 0.0;
            let p = (0..codebook.len()).filter_map(|feature_id| {
                let encoding = codebook.record(feature_id).codeword();
                if exact_count[(feature_id, 0)] < 100.0 || encoding[k] != bit {
                    None
                } else {
                    count += 1.0;
                    Some(feat_p[(feature_id, k)])
                }
            }).sum::<f64>() / count;
            if count == 0.0 {
                panic!("not enough counts to estimate error rate (at least 1 feature with 50 counts over all cells required)");
            }
            p
        };

        p0.push(Prob(p(&feat_p0, false)));
        p1.push(Prob(p(&feat_p1, true)));
    }

    (p0, p1)
}
