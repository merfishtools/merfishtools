#![allow(non_snake_case)]

use std::cmp;
use std::mem;
use std::collections::HashMap;
use std::cell::RefCell;

use rgsl::randist::multinomial::multinomial_pdf;
use itertools::Itertools;
use bit_vec::BitVec;

use bio::stats::combinatorics::combinations;
use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codebook, Codeword};
use ndarray::prelude::*;


/// Basic model for the probability of making i 1-0 and j 0-1 errors.
struct Xi {
    p0: Vec<LogProb>,
    p1: Vec<LogProb>,
    max_err: u8
}


impl Xi {
    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `p0` - error probabilities for 0-1 errors
    /// * `p1` - error probabilities for 1-0 errors
    /// * `max_err` - maximum number of errors considered
    fn new(p0: &[Prob], p1: &[Prob], max_err: u8) -> Self {
        let tolog = |&p| LogProb::from(p);
        Xi {
            p0: p0.iter().map(&tolog).collect_vec(),
            p1: p1.iter().map(&tolog).collect_vec(),
            max_err: max_err
        }
    }

    fn calc(&self, codeword: &Codeword, mask: Option<&BitVec>) -> Array2<LogProb> {
        let n = self.max_err as usize;
        let mut curr = Array::from_elem((n + 1, n + 1), LogProb::ln_zero());
        let mut prev = Array::from_elem((n + 1, n + 1), LogProb::ln_zero());
        prev[(1, 1)] = LogProb::ln_one();

        for k in 0..codeword.len() {
            for i in 1..curr.shape()[0] {
                for j in 1..curr.shape()[1] {

                    let (i_err, j_err, p_err) = if codeword.get(k).unwrap() {
                        (i - 1, j, self.p1[k])
                    } else {
                        (i, j - 1, self.p0[k])
                    };

                    let mut p = p_err.ln_one_minus_exp() + prev[(i, j)];

                    if mask.map_or(true, |mask| mask.get(k).unwrap()) {
                        p = p.ln_add_exp(p_err + prev[(i_err, j_err)]);
                    }

                    curr[(i, j)] = p;
                }
            }
            mem::swap(&mut prev, &mut curr);
        }

        prev.slice(s![1.., 1..]).to_owned()
    }
}


pub struct Params {
    codebook: Codebook,
    xi: Xi
}


#[derive(Debug)]
pub struct Events {
    exact: LogProb,
    mismatch: LogProb
}


impl Events {
    pub fn total(&self) -> LogProb {
        self.exact.ln_add_exp(self.mismatch)
    }
}


pub trait Model: Sync {

    fn params(&self) -> &Params;

    /// Probability to make exactly i 1->0 and j 0->1 errors towards any target.
    fn xi(&self, source: &codebook::Record, target: Option<&codebook::Record>) -> Array2<LogProb> {
        let mask = if let Some(target) = target {
            Some(source.diff(target))
        } else {
            None
        };

        self.params().xi.calc(source.codeword(), mask.as_ref())
    }

    /// Probability to see an exact readout or a readout with one mismatch given that we have a call.
    fn prob_call(&self, feature: &str) -> Events;

    /// Probability to see a miscalled exact readout or a miscalled readout with one mismatch.
    fn prob_miscall(&self, feature: &str) -> Events;

    // /// Probability to miss or miscall a readout because of too many errors.
    // fn prob_missed(&self, feature: &str) -> LogProb {
    //     self.prob_call_exact(feature).ln_add_exp(self.prob_call_mismatch(feature)).ln_one_minus_exp()
    // }
    //
    // fn prob_miscall(&self, feature: &str) -> LogProb {
    //     self.prob_miscall_exact(feature).ln_add_exp(self.prob_miscall_mismatch(feature))
    // }
    //
    // /// Probability to completely miss a readout (not miscalled).
    // fn prob_nocall(&self, feature: &str) -> LogProb {
    //     self.prob_miscall_exact(feature).ln_add_exp(self.prob_miscall_mismatch(feature)).ln_one_minus_exp()
    // }
}


pub struct MHD4 {
    params: Params
}


impl Model for MHD4 {

    fn params(&self) -> &Params {
        &self.params
    }

    fn prob_call(&self, feature: &str) -> Events {
        let xi = self.xi(self.params().codebook.get(feature), None);

        Events {
            exact: xi[(0, 0)],
            mismatch: xi[(1, 0)].ln_add_exp(xi[(0, 1)])
        }
    }

    fn prob_miscall(&self, feature: &str) -> Events {
        let feature = self.params().codebook.get(feature);

        let mut exact_summands = Vec::new();
        let mut mismatch_summands = Vec::new();
        for neighbor in self.params().codebook.neighbors(feature.name(), 4) {
            let xi = self.xi(feature, Some(&neighbor));
            exact_summands.push(xi[(2, 2)]);
            mismatch_summands.push(xi[(2, 1)]);
            mismatch_summands.push(xi[(1, 2)]);
            mismatch_summands.push(xi[(2, 3)]);
            mismatch_summands.push(xi[(3, 2)]);
        }

        Events {
            exact: LogProb::ln_sum_exp(&exact_summands),
            mismatch: LogProb::ln_sum_exp(&mismatch_summands)
        }
    }
}


pub struct MHD2 {
    params: Params
}


impl Model for MHD2 {

    fn params(&self) -> &Params {
        &self.params
    }

    fn prob_call(&self, feature: &str) -> Events {
        let xi = self.xi(self.params().codebook.get(feature), None);

        Events {
            exact: xi[(0, 0)],
            mismatch: LogProb::ln_zero()
        }
    }

    fn prob_miscall(&self, feature: &str) -> Events {
        let feature = self.params().codebook.get(feature);

        let mut exact_summands = Vec::new();
        for neighbor in self.params().codebook.neighbors(feature.name(), 4) {
            let xi = self.xi(feature, Some(&neighbor));
            exact_summands.push(xi[(2, 2)]);
        }
        for neighbor in self.params().codebook.neighbors(feature.name(), 4) {
            let xi = self.xi(feature, Some(&neighbor));
            exact_summands.push(xi[(1, 1)]);
        }

        Events {
            exact: LogProb::ln_sum_exp(&exact_summands),
            mismatch: LogProb::ln_zero()
        }
    }
}


pub fn new_model(p0: &[Prob], p1: &[Prob], mut codebook: Codebook) -> Box<Model> {
    // TODO get rid of hardcoded value for max_err
    let params = Params {codebook: codebook, xi: Xi::new(p0, p1, 6)};
    let model: Box<Model> = match params.codebook.min_dist {
        4 => Box::new(MHD4 { params: params }),
        2 => Box::new(MHD2 { params: params }),
        _ => panic!("Hamming distances other than 2 and 4 are unsupported at the moment.")
    };
    model
}


/// Readout probabilities.
pub struct Readout {
    prob_call_exact: Prob,
    prob_call_mismatch: Prob,
    prob_miscall_exact: Prob,
    prob_miscall_mismatch: Prob,
    prob_missed: Prob,
    prob_miscall: Prob,
    prob_call: Prob,
    margin: u32
}


impl Readout {
    pub fn new(feature: &str, mut model: &Box<Model>, window_width: u32) -> Self {
        let prob_call = model.prob_call(feature);
        let prob_miscall = model.prob_miscall(feature);

        Readout {
            prob_call_exact: Prob::from(prob_call.exact),
            prob_call_mismatch: Prob::from(prob_call.mismatch),
            prob_miscall_exact: Prob::from(prob_miscall.exact),
            prob_miscall_mismatch: Prob::from(prob_miscall.mismatch),
            prob_missed: Prob::from(prob_call.total().ln_one_minus_exp()),
            prob_miscall: Prob::from(prob_miscall.total()),
            prob_call: Prob::from(prob_miscall.total().ln_one_minus_exp()),
            margin: window_width / 2
        }
    }

    fn est_x(&self, n: f64) -> u32 {
        (
            n * *(self.prob_call_exact * self.prob_call) +
            n * *(self.prob_call_mismatch * self.prob_call) +
            n * *(self.prob_missed * self.prob_call)
        ).round() as u32
    }

    /// Upper and lower bound for MAP
    pub fn map_bounds(&self, count: u32, count_exact: u32) -> (u32, u32) {
        // estimate n from exact readouts
        let n_0 = count_exact as f64 / *(
            self.prob_call_exact * self.prob_call + self.prob_miscall_exact
        );
        let n_1 = {
            // estimate n from corrected readouts
            let p = self.prob_call_mismatch * self.prob_call + self.prob_miscall_mismatch;
            if *p > 0.0 {
                let n_1 = (count - count_exact) as f64 / *p;
                n_1
            } else {
                n_0
            }
        };

        // estimate lower and upper bound of x
        // with MHD2, this is the same and the estimate is equivalent to the MAP
        let mut x_0 = self.est_x(n_0);
        let mut x_1 = self.est_x(n_1);
        if x_0 > x_1 {
            mem::swap(&mut x_0, &mut x_1);
        }

        (x_0 as u32, x_1 as u32)
    }

    pub fn naive_estimate(&self, count: u32) -> u32 {
        let n_avg = count as f64 / *(
            self.prob_call_exact * self.prob_call +
            self.prob_call_mismatch * self.prob_call +
            self.prob_miscall_exact +
            self.prob_miscall_mismatch
        );
        self.est_x(n_avg) as u32
    }

    /// Prior window for calculating expression PMF
    pub fn window(&self, count: u32, count_exact: u32) -> (u32, u32) {
        let x = self.naive_estimate(count);

        (cmp::max(x as i32 - self.margin as i32, 0) as u32, x + self.margin)
    }

    pub fn likelihood(&self, x: u32, count: u32, count_exact: u32) -> LogProb {
        let count = count;
        let count_exact = count_exact;
        assert!(count >= count_exact);

        let mut summands = Vec::new();
        let probs = [
            *(self.prob_call_exact * self.prob_call), // Pr(H=0, E=e) = Pr(H=0 | E=e) * Pr(E=e)
            *(self.prob_call_mismatch * self.prob_call),
            *(self.prob_missed * self.prob_call),
            *self.prob_miscall_exact,
            *self.prob_miscall_mismatch
        ];

        for i in 0..(cmp::min(x, count) + 1) {
            let jmax = cmp::min(count_exact, i);
            // i - j <= count - count_exact
            let jmin = if count_exact + i > count { count_exact + i - count } else { 0 };
            for j in jmin..(jmax + 1) {
                // let n = combinations((count - count_exact) as u64, (i - j) as u64).ln() + combinations(count_exact as u64, j as u64).ln();
                let k = x - i;
                let exact_miscalls = count_exact - j;
                let mismatch_miscalls = count - count_exact - (i - j);

                let p = LogProb::from(Prob(
                    multinomial_pdf(&probs, &[j, i - j, k, exact_miscalls, mismatch_miscalls])
                ));
                summands.push(p);
            }
        }
        let likelihood = LogProb::ln_sum_exp(&summands);
        assert!(!likelihood.is_nan());
        likelihood
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use super::*;
    use nalgebra::ApproxEq;
    use io;
    use bio::stats::combinatorics::combinations;

    fn setup_mhd4() -> Box<Model> {
        new_model(
            &[Prob(0.04); 16],
            &[Prob(0.1); 16],
            io::codebook::Codebook::from_file("tests/codebook/simulated-MHD4.txt").unwrap()
        )
    }


    fn setup_mhd2() -> Box<Model> {
        new_model(
            &[Prob(0.04); 16],
            &[Prob(0.1); 16],
            io::codebook::Codebook::from_file("tests/codebook/simulated-MHD2.txt").unwrap()
        )
    }


    #[test]
    fn test_prob_call_exact() {
        let feat = "COL5A1";
        let model = setup_mhd4();
        let p = model.prob_call(feat).exact;
        println!("{}", *p);
        assert!(p.exp().approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_prob_call_mismatch() {
        let feat = "COL5A1";
        let model = setup_mhd4();
        let p = model.prob_call(feat).mismatch;
        println!("{}", *p);
        assert!(p.exp().approx_eq(&0.3796656011293904));
    }

    #[test]
    fn test_prob_miscall_exact() {
        let feat = "COL5A1";
        let model = setup_mhd4();
        let p = model.prob_miscall(feat).exact;
        assert_relative_eq!(p.exp(), 0.00031018431464819473)
    }

    #[test]
    fn test_prob_miscall_mismatch() {
        let feat = "COL5A1";
        let model = setup_mhd4();
        let p = model.prob_miscall(feat).mismatch;
        println!("{}", *p);
        assert_relative_eq!(p.exp(), 0.0204, epsilon=0.001);
    }

    #[test]
    fn test_prob_missed() {
        let feat = "COL5A1";
        let model = setup_mhd4();
        let p = model.prob_call(feat).total().ln_one_minus_exp();
        assert_relative_eq!(p.exp(), 0.21833552708654924);
    }

    #[test]
    fn test_mhd2() {
        let feat = "COL7A1";
        let model = setup_mhd2();
        println!("{:?}", model.prob_call(feat));
        println!("{:?}", model.prob_miscall(feat));
    }

    fn comb(i: u8, j: u8) -> f64 {
        combinations(4, i as u64) * combinations(12, j as u64)
    }

    #[test]
    fn test_xi00() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(0, 0)].exp();
        assert_relative_eq!(p, comb(0,0) * 0.4019988717840602);
    }

    #[test]
    fn test_xi10() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(1, 0)].exp();
        assert_relative_eq!(p, comb(1,0) * 0.04466654130934002);
    }

    #[test]
    fn test_xi01() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(0, 1)].exp();
        assert_relative_eq!(p, comb(0, 1) * 0.01674995299100251);
    }

    #[test]
    fn test_xi22() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(2, 2)].exp();
        assert_relative_eq!(p, comb(2, 2) * 8.616230962449852e-06);
    }

    #[test]
    fn test_xi21() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(2, 1)].exp();
        assert_relative_eq!(p, comb(2, 1) * 0.00020678954309879646);
    }

    #[test]
    fn test_xi12() {
        let model = setup_mhd4();
        let feat = model.params().codebook.get("COL5A1");
        let p = model.xi(feat, None)[(1, 2)].exp();
        assert_relative_eq!(p, comb(1, 2) * 7.754607866204867e-05);
    }

    #[test]
    fn test_window() {
        let model = setup_mhd4();
        let readout = Readout::new("COL5A1", &model, 100);
        let (lower, upper) = readout.window(175, 75);
        println!("{} {}", lower, upper);
        for x in lower..upper {
            println!("{}={}", x, *readout.likelihood(x, 175, 75));
        }
        assert!(readout.likelihood(lower, 175, 75).exp().approx_eq(&0.0));
        assert!(readout.likelihood(upper, 175, 75).exp().approx_eq(&0.0));
    }
}
