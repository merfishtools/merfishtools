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


pub fn estimate<I: Iterator<Item=(String, Readout)>>(
    codebook: &Codebook,
    readouts: I
) -> (Vec<Prob>, Vec<Prob>) {
    let max_diff = 0.001;
    let shape = (codebook.len(), codebook.N as usize);

    let mut observed_01_errors = Array2::from_elem(shape, 0.0);
    let mut observed_10_errors = Array2::from_elem(shape, 0.0);
    let mut feat_count = Array2::from_elem((codebook.len(), 1), 0.0);
    for (feat, readout) in readouts {
        let feature_id = codebook.get_id(&feat);
        let encoding = codebook.record(feature_id).codeword();

        // count features
        feat_count[(feature_id, 0)] += 1.0;

        for k in 0..codebook.N as usize {
            if encoding[k] == false {
                // count 01 errors
                if readout[k] == true {
                    observed_01_errors[(feature_id, k)] += 1.0;
                }
            } else {
                // count 10 errors
                if readout[k] == false {
                    observed_10_errors[(feature_id, k)] += 1.0;
                }
            }
        }
    }

    // calculate observed (joint) error rates
    let observed_p0 = observed_01_errors / &feat_count;
    let observed_p1 = observed_10_errors / &feat_count;

    // run EM algorithm to estimate error rates
    let mut prob_no_error = Array2::from_elem(
        (codebook.len(), codebook.N as usize), LogProb::ln_one()
    );


    let mut p0 = Array1::random(codebook.N as usize, Range::new(0.0, 1.0));
    let mut p1 = Array1::random(codebook.N as usize, Range::new(0.0, 1.0));

    let mut p0_last = Array1::from_elem(codebook.N as usize, 0.0);
    let mut p1_last = Array1::from_elem(codebook.N as usize, 0.0);
    for i in 0..1000 {
        // E-step: estimate probability to have no errors in rest of readout
        for feature_id in 0..codebook.len() {
            let encoding = codebook.record(feature_id).codeword();
            for k in 0..codebook.N as usize {
                let p = LogProb((0..codebook.N as usize).filter_map(|l| {
                    if l != k {
                        let p_err = if encoding[l] == false {
                            p0[k]
                        } else {
                            p1[k]
                        };
                        Some(*LogProb((p_err as f64).ln()).ln_one_minus_exp())
                    } else { None }
                }).sum());
                prob_no_error[(feature_id, k)] = p;
            }
        }

        // M-step: find p0 and p1 that maximizes the likelihood
        mem::swap(&mut p0, &mut p0_last);
        mem::swap(&mut p1, &mut p1_last);
        for k in 0..codebook.N as usize {
            // solve least squares: take mean of all estimates
            let p = |observed_p: &Array2<f64>, bit: bool| {
                let mut count = 0.0;
                let p = (0..codebook.len()).filter_map(|feature_id| {
                    let encoding = codebook.record(feature_id).codeword();
                    if feat_count[(feature_id, 0)] < 10.0 || encoding[k] != bit {
                        None
                    } else {
                        count += 1.0;
                        let p = observed_p[(feature_id, k)] * prob_no_error[(feature_id, k)].exp();
                        Some(p)
                    }
                }).sum::<f64>() / count;
                if count == 0.0 {
                    panic!("not enough counts to estimate error rate");
                }
                p
            };
            //debug!("----------- p0 ------------");
            p0[k] = p(&observed_p0, false);
            //debug!("DEBUG: p0={}", p0[k]);
            //debug!("----------- p1 ------------");
            p1[k] = p(&observed_p1, true);
            //debug!("DEBUG: p1={}", p1[k]);
        }

        // check convergence
        let abs = |d: &f64| NotNaN::new(d.abs()).unwrap();
        //debug!("{:?}", p0);
        let convergence = p0.all_close(&p0_last, max_diff) && p1.all_close(&p1_last, max_diff);

        if convergence && i >= 5 {
            debug!("convergence reached");
            debug!("{:?}", p1);
            break;
        }

        debug!("EM-algorithm step {} of {}.", i + 1, 100);
    }

    (
        p0.iter().map(|p| Prob(*p)).collect_vec(),
        p1.iter().map(|p| Prob(*p)).collect_vec()
    )
}
