#![allow(non_snake_case)]

use std::cmp;
use std::mem;

use rgsl::randist::multinomial::multinomial_pdf;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs::{Prob, LogProb};
use bio::stats::logprobs;

use io::codebook::Codebook;

pub struct Params {
    N: u8,
    m: u8,
    p0: Prob,
    p1: Prob,
    codebook: Codebook
}


pub trait Model: Sync {

    fn params(&self) -> &Params;

    /// Probability to make exactly i 1->0 and j 0->1 errors.
    fn xi(&self, i: u8, j: u8) -> Prob {
        self.params().p1.powi(i as i32) * self.params().p0.powi(j as i32) *
        (1.0 - self.params().p1).powi((self.params().m - i) as i32) * (1.0 - self.params().p0).powi((self.params().N - self.params().m - j) as i32)
    }

    /// Number of possibilities to select i ones and j zeros.
    fn psi(&self, i: u8, j: u8) -> f64 {
        combinations(self.params().m as u64, i as u64) * combinations((self.params().N - self.params().m) as u64, j as u64)
    }

    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self) -> Prob;

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self) -> Prob;

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> Prob;

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self, feature: &str) -> Prob;

    /// Probability to miss or miscall a readout because of too many errors.
    fn prob_missed(&self) -> Prob {
        1.0 - self.prob_call_exact() - self.prob_call_mismatch()
    }

    fn prob_miscall(&self, feature: &str) -> Prob {
        self.prob_miscall_exact(feature) + self.prob_miscall_mismatch(feature)
    }

    /// Probability to completely miss a readout (not miscalled).
    fn prob_nocall(&self, feature: &str) -> Prob {
        1.0 - self.prob_miscall_exact(feature) - self.prob_miscall_mismatch(feature)
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
    fn prob_call_exact(&self) -> Prob {
        self.xi(0, 0)
    }

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self) -> Prob {
        (self.params().m as f64 * self.xi(1, 0) + (self.params().N - self.params().m) as f64 * self.xi(0, 1))
    }

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> Prob {
        self.params().codebook.neighbors(feature, 4).len() as f64 * self.xi(2, 2)
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self, feature: &str) -> Prob {
        let n = self.params().codebook.neighbors(feature, 4).len() as f64;

        n * 2 as f64 * self.xi(2, 1) + n * 2 as f64 * self.xi(1, 2) +
        n * (self.params().m - 2) as f64 * self.xi(3, 2) + n * (self.params().N - self.params().m - 2) as f64 * self.xi(2, 3)
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
    fn prob_call_exact(&self) -> Prob {
        self.xi(0, 0)
    }

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self) -> Prob {
        0.0
    }

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self, feature: &str) -> Prob {
        self.params().codebook.neighbors(feature, 2).len() as f64 * self.xi(1, 1) +
        self.params().codebook.neighbors(feature, 4).len() as f64 * self.xi(2, 2)
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    #[allow(unused_variables)]
    fn prob_miscall_mismatch(&self, feature: &str) -> Prob {
        0.0
    }
}


pub fn new_model(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8, codebook: Codebook) -> Box<Model> {
    let params = Params { N: N, m: m, p0: p0, p1: p1, codebook: codebook};
    let model: Box<Model> = match dist {
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
            prob_call_exact: model.prob_call_exact(),
            prob_call_mismatch: model.prob_call_mismatch(),
            prob_miscall_exact: model.prob_miscall_exact(feature),
            prob_miscall_mismatch: model.prob_miscall_mismatch(feature),
            prob_missed: model.prob_missed(),
            prob_nocall: model.prob_nocall(feature),
            prob_miscall: prob_miscall,
            prob_call: 1.0 - prob_miscall,
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
                //let n = combinations((count - count_exact) as u64, (i - j) as u64).ln() + combinations(count_exact as u64, j as u64).ln();
                let k = x - i;
                let exact_miscalls = count_exact - j;
                let mismatch_miscalls = count - count_exact - (i - j);

                let p = multinomial_pdf(&probs, &[j, i - j, k, exact_miscalls, mismatch_miscalls]).ln();
                summands.push(p);
            }
        }
        let likelihood = logprobs::sum(&summands);
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
        new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("test/codebook/simulated-MHD4.txt", 4).unwrap().codebook())
    }


    fn setup_mhd2() -> Box<Model> {
        new_model(14, 4, 0.04, 0.1, 2, io::codebook::Reader::from_file("test/codebook/simulated-MHD2.txt", 2).unwrap().codebook())
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
