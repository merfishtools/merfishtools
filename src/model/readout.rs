#![allow(non_snake_case)]

use std::cmp;

use rgsl::randist::multinomial::multinomial_pdf;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs::{Prob, LogProb};
use bio::stats::logprobs;


struct Params {
    N: u8,
    m: u8,
    p0: Prob,
    p1: Prob
}


trait Model {

    fn params(&self) -> &Params;

    /// Probability to make exactly i 1->0 and j 0->1 errors.
    fn xi(&self, i: u8, j: u8) -> Prob {
        self.params().p1.powi(i as i32) * self.params().p0.powi(j as i32) *
        (1.0 - self.params().p1).powi((self.params().m - i) as i32) * (1.0 - self.params().p0).powi((self.params().N - self.params().m - j) as i32)
    }

    /// Number of possibilities for i 1->0 and j 0->1 errors.
    fn psi(&self, i: u8, j: u8) -> f64 {
        combinations(self.params().m as u64, i as u64) * combinations((self.params().N - self.params().m) as u64, j as u64)
    }

    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self) -> Prob;

    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self) -> Prob;

    /// Probability to see an exact readout given that we have a miscall.
    fn prob_miscall_exact(&self) -> Prob;

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self) -> Prob;

    /// Probability to miss or miscall a readout because of too many errors.
    fn prob_missed(&self) -> Prob {
        1.0 - self.prob_call_exact() - self.prob_call_mismatch()
    }

    fn prob_miscall(&self) -> Prob {
        self.prob_miscall_exact() + self.prob_miscall_mismatch()
    }

    /// Probability to completely miss a readout (not miscalled).
    fn prob_nocall(&self) -> Prob {
        1.0 - self.prob_miscall_exact() - self.prob_miscall_mismatch()
    }
}


struct MHD4 {
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
    fn prob_miscall_exact(&self) -> Prob {
        self.psi(2,2) * self.xi(2, 2)
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self) -> Prob {
        self.psi(2, 1) * self.xi(2, 1) + self.psi(1, 2) * self.xi(1, 2) +
        self.psi(2, 3) * self.xi(2, 3) + self.psi(3, 2) * self.xi(3, 2)
    }
}


struct MHD2 {
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
    fn prob_miscall_exact(&self) -> Prob {
        (1..4).fold(0.0, |p, i| p + self.psi(i, i) * self.xi(i, i))
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self) -> Prob {
        0.0
    }
}


/// Readout probabilities.
pub struct Readout {
    prob_call_exact: Prob,
    prob_call_mismatch: Prob,
    prob_miscall_exact: Prob,
    prob_miscall_mismatch: Prob,
    prob_missed: Prob,
    prob_nocall: Prob,
    prob_miscall: Prob
}


impl Readout {
    pub fn new(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8) -> Self {
        let params = Params { N: N, m: m, p0: p0, p1: p1};
        let model: Box<Model> = match dist {
            4 => Box::new(MHD4 { params: params }),
            2 => Box::new(MHD2 { params: params }),
            _ => panic!("Hamming distances other than 2 and 4 are unsupported at the moment.")
        };
        Readout {
            prob_call_exact: model.prob_call_exact(),
            prob_call_mismatch: model.prob_call_mismatch(),
            prob_miscall_exact: model.prob_miscall_exact(),
            prob_miscall_mismatch: model.prob_miscall_mismatch(),
            prob_missed: model.prob_missed(),
            prob_nocall: model.prob_nocall(),
            prob_miscall: model.prob_miscall()
        }
    }

    pub fn window(&self, count: u32) -> (u32, u32) {
        let prob_call = 1.0 - self.prob_miscall;
        //let n_0 = (count - count_exact) as f64 / (self.prob_call_mismatch * (1.0 - self.prob_miscall) + self.prob_miscall_mismatch);
        //let n_1 = count_exact as f64 / (self.prob_call_exact * (1.0 - self.prob_miscall) + self.prob_miscall_exact);
        let n = count as f64 / (self.prob_call_exact * prob_call +  self.prob_call_mismatch * prob_call + self.prob_miscall_exact + self.prob_miscall_mismatch);
        let x = (n * self.prob_call_exact * prob_call + n * self.prob_call_mismatch * prob_call + n * self.prob_missed).round() as i32;

        (cmp::max(x - 50, 0) as u32, x as u32 + 50)
    }

    pub fn likelihood(&self, x: u32, count: u32, count_exact: u32) -> LogProb {
        let x = x;
        let count = count;
        let count_exact = count_exact;
        assert!(count >= count_exact);

        let mut summands = Vec::new();
        let prob_call = 1.0 - self.prob_miscall; // Pr(E=e) = 1 - Pr(E!=e)
        let probs = [
            self.prob_call_exact * prob_call, // Pr(H=0, E=e) = Pr(H=0 | E=e) * Pr(E=e)
            self.prob_call_mismatch * prob_call,
            self.prob_missed * prob_call,
            self.prob_miscall_exact,
            self.prob_miscall_mismatch
        ];

        for i in 0..(cmp::min(x, count) + 1) {
            let jmax = cmp::min(count_exact, i);
            let jmin = if count_exact + i > count { count_exact + i - count } else { 0 };
            for j in jmin..(jmax + 1) {
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

    use super::{MHD4, MHD2, Params, Model, Readout};
    use nalgebra::ApproxEq;


    const factory: MHD4 = MHD4 { params: Params { N: 16, m: 4, p0: 0.04, p1: 0.1 } };


    #[test]
    fn test_prob_call_exact() {
        let p = factory.prob_call_exact();
        println!("{}", p);
        assert!(p.approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_prob_call_mismatch() {
        let p = factory.prob_call_mismatch();
        println!("{}", p);
        assert!(p.approx_eq(&0.3796656011293904));
    }

    #[test]
    fn test_prob_miscall_exact() {
        let p = factory.prob_miscall_exact();
        println!("{}", p);
        assert!(p.approx_eq(&0.003412027461130142));
    }

    #[test]
    fn test_prob_miscall_mismatch() {
        let p = factory.prob_miscall_mismatch();
        println!("{}", p);
        assert!(p.approx_eq(&0.03608764734772748));
    }

    #[test]
    fn test_prob_missed() {
        let p = factory.prob_missed();
        println!("{}", p);
        assert!(p.approx_eq(&0.21833552708654924));
    }

    #[test]
    fn test_mhd2() {
        let model = MHD2 { params: Params { N: 14, m: 4, p0: 0.04, p1: 0.1 } };
        println!("{}", model.prob_call_exact());
        println!("{}", model.prob_call_mismatch());
        println!("{}", model.prob_miscall_exact());
        println!("{}", model.prob_miscall_mismatch());
        println!("{}", model.prob_missed());
    }

    #[test]
    fn test_xi() {
        let p = factory.xi(0, 0);
        assert!(p.approx_eq(&0.4019988717840602));
        let p = factory.xi(1, 0);
        assert!(p.approx_eq(&0.04466654130934002));
        let p = factory.xi(0, 1);
        assert!(p.approx_eq(&0.01674995299100251));
        let p = factory.xi(2, 2);
        assert!(p.approx_eq(&8.616230962449852e-06));
        let p = factory.xi(2, 1);
        assert!(p.approx_eq(&0.00020678954309879646));
        let p = factory.xi(1, 2);
        assert!(p.approx_eq(&7.754607866204867e-05));
    }

    #[test]
    fn test_window() {
        let readout = Readout::new(16, 4, 0.04, 0.1, 4);
        let (lower, upper) = readout.window(175);
        println!("{} {}", lower, upper);
        assert!(readout.likelihood(lower, 175, 25).exp().approx_eq(&0.0));
        assert!(readout.likelihood(upper, 175, 25).exp().approx_eq(&0.0));
    }
}
