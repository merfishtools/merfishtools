use std::f64;

use ndarray::prelude::*;

use bio::stats::Prob;

use io::codebook::Codebook;
use io::merfishdata::Readout;


/// Estimate position-wise error rates as presented by Chen et al. 2015.
pub fn estimate<I: Iterator<Item = (String, Readout)>>(
    codebook: &Codebook,
    readouts: I,
) -> (Vec<Prob>, Vec<Prob>) {
    let _max_diff = 0.001;
    let shape = (codebook.len(), codebook.N as usize);

    let mut observed_01_errors = Array2::from_elem(shape, 0.0);
    let mut observed_10_errors = Array2::from_elem(shape, 0.0);
    let mut exact_count = Array2::from_elem((codebook.len(), 1), 0.0);
    let mut feat_count = Array1::from_elem(codebook.len(), 0);

    for (feat, readout) in readouts {
        let feature_id = codebook.get_id(&feat);
        let rec = codebook.record(feature_id);
        if !rec.expressed() {
            continue;
        }
        let encoding = rec.codeword();
        feat_count[feature_id] += 1;

        let mut err_count = 0;
        for k in 0..codebook.N as usize {
            match (encoding[k], readout[k]) {
                (false, true) => {
                    observed_01_errors[(feature_id, k)] += 1.0;
                    err_count += 1;
                }
                (true, false) => {
                    observed_10_errors[(feature_id, k)] += 1.0;
                    err_count += 1;
                }
                _ => (),
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
            let mut total_weight = 0.0;
            let p = (0..codebook.len())
                .filter_map(|feature_id| {
                    let encoding = codebook.record(feature_id).codeword();
                    if encoding[k] != bit || exact_count[(feature_id, 0)] == 0.0 {
                        None
                    } else {
                        let w = feat_count[feature_id] as f64;
                        total_weight += w;
                        Some(feat_p[(feature_id, k)] * w)
                    }
                })
                .sum::<f64>() / total_weight;
            if total_weight == 0.0 {
                panic!("no readouts to estimate error rate at position {}", k);
            }
            p
        };

        p0.push(Prob(p(&feat_p0, false)));
        p1.push(Prob(p(&feat_p1, true)));
    }

    (p0, p1)
}
