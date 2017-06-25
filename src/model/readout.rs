#![allow(non_snake_case)]

use std::cmp;
use std::mem;
use std::collections::HashMap;

use rgsl::randist::multinomial::multinomial_pdf;

use bio::stats::combinatorics::combinations;
use bio::stats::{Prob, LogProb};

use io::codebook::{Codebook, Codeword};
use ndarray::prelude::*;


/// Basic model for the probability of making i 1-0 and j 0-1 errors.
struct Xi {
    p0: Vec<LogProb>,
    p1: Vec<LogProb>,
    curr: Array2<LogProb>,
    prev: Array2<LogProb>
}


impl Xi {
    fn new(p0: &[Prob], p1: &[Prob], max_err: u8) -> Self {
        let tolog = |p| LogProb::from(p);
        let n = max_err as usize;
        Xi {
            curr: Array::default((n + 1, n + 1)),
            prev: Array::default((n + 1, n + 1)),
            p0: p0.iter().map(tolog).collect_vec(),
            p1: p1.iter.map(tolog).collect_vec()
        }
    }

    fn calc(&mut self, codeword: &Codeword) -> Array2<LogProb> {
        self.curr[(.., ..)] = LogProb::ln_zero();
        self.prev[(.., ..)] = LogProb::ln_zero();

        self.prev[(0, 0)] = LogProb::ln_one();

        for k in 0..codeword.len() {
            for i in 1..self.curr.shape()[0] {
                for j in 1..self.curr.shape()[1] {
                    self.curr[(i, j)] = if codeword.get(k) {
                        (self.p1[k] + self.prev[(i - 1, j)]).ln_add_exp(
                            self.p1[k].ln_one_minus_exp() + self.prev[(i, j)]
                        )
                    } else {
                        (self.p0[k] + self.prev[(i, j - 1)]).ln_add_exp(
                            self.p0[k].ln_one_minus_exp() + self.prev[(i, j)]
                        )
                    };
                }
            }
            mem::swap(&mut self.prev, &mut self.curr);
        }

        self.prev.clone()
    }
}


pub struct Params {
    codebook: Codebook,
    xi: HashMap<String, Array2<LogProb>>
}


pub trait Model: Sync {

    fn params(&self) -> &Params;

    /// Probability to make exactly i 1->0 and j 0->1 errors.
    fn xi(&self, feature: &str, i: i8, j: i8) -> LogProb {
        self.xi.get(feature).unwrap()[(i, j)]
    }

    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self, feature: &str) -> LogProb;

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self, feature: &str) -> LogProb;

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> LogProb;

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self, feature: &str) -> LogProb;

    /// Probability to miss or miscall a readout because of too many errors.
    fn prob_missed(&self, feature: &str) -> LogProb {
        self.prob_call_exact(feature).ln_add_exp(self.prob_call_mismatch(feature)).ln_one_minus_exp()
    }

    fn prob_miscall(&self, feature: &str) -> LogProb {
        self.prob_miscall_exact(feature).ln_add_exp(self.prob_miscall_mismatch(feature))
    }

    /// Probability to completely miss a readout (not miscalled).
    fn prob_nocall(&self, feature: &str) -> LogProb {
        self.prob_miscall_exact(feature).ln_add_exp(self.prob_miscall_mismatch(feature)).ln_one_minus_exp()
    }
}


pub struct MHD4 {
    params: Params
}


impl Model for MHD4 {

    fn params(&self) -> &Params {
        &self.params
    }

    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self, feature: &str) -> LogProb {
        self.xi(feature, 0, 0)
    }

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self, feature: &str) -> LogProb {
        self.xi(feature, 1, 0).ln_add_exp(self.xi(feature, 0, 1))
    }

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> LogProb {
        LogProb::ln_sum_exp(&self.params().codebook.neighbors(feature, 4).map(|neighbor| {
            self.xi(neighbor.name(), 2, 2)
        }).collect_vec())
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self, feature: &str) -> Prob {
        let neighbors = self.params().codebook.neighbors(feature, 4);
        let mut summands = Vec::with_capacity(neighbors.len() * 4);
        for neighbor in neighbors {
            summands.push(self.xi(neighbor.name(), 2, 1));
            summands.push(self.xi(neighbor.name(), 1, 2));
            summands.push(self.xi(neighbor.name(), 2, 3));
            summands.push(self.xi(neighbor.name(), 3, 2));
        }

        LogProb::ln_sum_exp(&summands)
    }
}


pub struct MHD2 {
    params: Params
}


impl Model for MHD2 {

    fn params(&self) -> &Params {
        &self.params
    }

    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self, feature: &str) -> LogProb {
        self.xi(feature, 0, 0)
    }

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self, feature: &str) -> LogProb {
        LogProb::ln_zero()
    }

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> Prob {
        let neighbors_4 = self.params().codebook.neighbors(feature, 4);
        let neighbors_2 = self.params().codebook.neighbors(feature, 2);
        let mut summands = Vec::with_capacity(neighbors_4.len() + neighbors_2.len());
        for neighbor in neighbors_4 {
            summands.push(self.xi(neighbor.name(), 2, 2));
        }
        for neighbor in neighbors_2 {
            summands.push(self.xi(neighbor.name(), 1, 1));
        }

        LogProb::ln_sum_exp(&summands)
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    #[allow(unused_variables)]
    fn prob_miscall_mismatch(&self, feature: &str) -> LogProb {
        LogProb::ln_zero()
    }
}


pub fn new_model(p0: &[Prob], p1: &[Prob], codebook: Codebook) -> Box<Model> {
    let mut xi = HashMap::new();
    let mut xicalc = Xi::new(p0, p1);
    for rec in codebook.records() {
        xi.insert(rec.name, xicalc.calc(rec.codeword));
    }

    let params = Params {codebook: codebook, xi: xi};
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
    prob_nocall: Prob,
    prob_miscall: Prob,
    prob_call: Prob,
    margin: u32
}


impl Readout {
    pub fn new(feature: &str, model: &Box<Model>, window_width: u32) -> Self {
        let prob_miscall = model.prob_miscall(feature);
        Readout {
            prob_call_exact: Prob::from(model.prob_call_exact()),
            prob_call_mismatch: Prob::from(model.prob_call_mismatch()),
            prob_miscall_exact: Prob::from(model.prob_miscall_exact(feature)),
            prob_miscall_mismatch: Prob::from(model.prob_miscall_mismatch(feature)),
            prob_missed: Prob::from(model.prob_missed()),
            prob_nocall: Prob::from(model.prob_nocall(feature)),
            prob_miscall: Prob::from(prob_miscall),
            prob_call: Prob::from(prob_miscall.ln_one_minus_exp()),
            margin: window_width / 2
        }
    }

    fn est_x(&self, n: f64) -> u32 {
        (
            n * self.prob_call_exact * self.prob_call +
            n * self.prob_call_mismatch * self.prob_call +
            n * self.prob_missed * self.prob_call
        ).round() as u32
    }

    /// Upper and lower bound for MAP
    pub fn map_bounds(&self, count: u32, count_exact: u32) -> (u32, u32) {
        // estimate n from exact readouts
        let n_0 = count_exact as f64 / (self.prob_call_exact * self.prob_call + self.prob_miscall_exact);
        let n_1 = {
            // estimate n from corrected readouts
            let p = self.prob_call_mismatch * self.prob_call + self.prob_miscall_mismatch;
            if p > 0.0 {
                let n_1 = (count - count_exact) as f64 / p;
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
        let n_avg = count as f64 / (
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
            self.prob_call_exact * self.prob_call, // Pr(H=0, E=e) = Pr(H=0 | E=e) * Pr(E=e)
            self.prob_call_mismatch * self.prob_call,
            self.prob_missed * self.prob_call,
            self.prob_miscall_exact,
            self.prob_miscall_mismatch
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

    pub fn expected_total(&self, total_count: u32) -> f64 {
        // c = x - x * p_nocall
        // x = c / (1 - p_nocall)
        total_count as f64 / (1.0 - self.prob_nocall)
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use super::*;
    use nalgebra::ApproxEq;
    use io;

    fn setup_mhd4() -> Box<Model> {
        new_model(0.04, 0.1, io::codebook::Codebook::from_file("test/codebook/simulated-MHD4.txt").unwrap())
    }


    fn setup_mhd2() -> Box<Model> {
        new_model(0.04, 0.1, io::codebook::Codebook::from_file("test/codebook/simulated-MHD2.txt").unwrap())
    }


    #[test]
    fn test_prob_call_exact() {
        let model = setup_mhd4();
        let p = model.prob_call_exact();
        println!("{}", p);
        assert!(p.approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_prob_call_mismatch() {
        let model = setup_mhd4();
        let p = model.prob_call_mismatch();
        println!("{}", p);
        assert!(p.approx_eq(&0.3796656011293904));
    }

    #[test]
    fn test_prob_miscall_exact() {
        let model = setup_mhd4();
        let p = model.prob_miscall_exact("COL5A1");
        println!("{}", p);
        assert!(p.approx_eq(&0.00031018431464819473));
    }

    #[test]
    fn test_prob_miscall_mismatch() {
        let model = setup_mhd4();
        let p = model.prob_miscall_mismatch("COL5A1");
        println!("{}", p);
        assert!(p.approx_eq(&0.020670338078917206));
    }

    #[test]
    fn test_prob_missed() {
        let model = setup_mhd4();
        let p = model.prob_missed();
        println!("{}", p);
        assert!(p.approx_eq(&0.21833552708654924));
    }

    #[test]
    fn test_mhd2() {
        let model = setup_mhd2();
        println!("{}", model.prob_call_exact());
        println!("{}", model.prob_call_mismatch());
        println!("{}", model.prob_miscall_exact("COL7A1"));
        println!("{}", model.prob_miscall_mismatch("COL7A1"));
        println!("{}", model.prob_missed());
    }

    #[test]
    fn test_xi() {
        let model = setup_mhd4();
        let p = model.xi(0, 0);
        assert!(p.approx_eq(&0.4019988717840602));
        let p = model.xi(1, 0);
        assert!(p.approx_eq(&0.04466654130934002));
        let p = model.xi(0, 1);
        assert!(p.approx_eq(&0.01674995299100251));
        let p = model.xi(2, 2);
        assert!(p.approx_eq(&8.616230962449852e-06));
        let p = model.xi(2, 1);
        assert!(p.approx_eq(&0.00020678954309879646));
        let p = model.xi(1, 2);
        assert!(p.approx_eq(&7.754607866204867e-05));
    }

    #[test]
    fn test_window() {
        let model = setup_mhd4();
        let readout = Readout::new("COL5A1", &model, 100);
        let (lower, upper) = readout.window(175, 75);
        println!("{} {}", lower, upper);
        for x in lower..upper {
            println!("{}={}", x, readout.likelihood(x, 175, 75));
        }
        assert!(readout.likelihood(lower, 175, 75).exp().approx_eq(&0.0));
        assert!(readout.likelihood(upper, 175, 75).exp().approx_eq(&0.0));
    }
}
