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

    /// Probability to completely miss a readout because of too many errors.
    fn prob_missed(&self) -> Prob;
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

    /// Probability to completely miss a readout because of too many errors.
    fn prob_missed(&self) -> Prob {
        let mut p = 0.0;
        for k in 2..6 {
            for i in 0..cmp::min(self.params().m, k) + 1 {
                p += self.psi(i, k - i) * self.xi(i, k - i);
            }
        }
        p
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
        (1..3).fold(0.0, |p, i| p + self.psi(i, i) * self.xi(i, i))
    }

    /// Probability to see a readout with one mismatch given that we have a miscall.
    fn prob_miscall_mismatch(&self) -> Prob {
        0.0
    }

    /// Probability to miss a readout because of too many errors.
    fn prob_missed(&self) -> Prob {
        let mut p = 0.0;
        for k in 1..6 {
            for i in 0..cmp::min(self.params().m, k) + 1 {
                p += self.psi(i, k - i) * self.xi(i, k - i);
            }
        }
        p
    }
}


/// Readout probabilities.
pub struct Readout {
    prob_call_exact: Prob,
    prob_call_mismatch: Prob,
    prob_miscall_exact: LogProb,
    prob_miscall_mismatch: LogProb,
    prob_missed: Prob
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
            prob_miscall_exact: model.prob_miscall_exact().ln(),
            prob_miscall_mismatch: model.prob_miscall_mismatch().ln(),
            prob_missed: model.prob_missed()
        }
    }

    pub fn window(&self, count: u32) -> (u32, u32) {
        // i = count and j = count_exact maximize the inner probability of the likelihood because
        // we have no miscalls then.
        // Hence we can estimate a good center using the expected value of the multinomial distribution.
        // E(j) = x * Pr(H=0 | E=e)
        // E(i - j) = x * Pr(H=1 | E=e)
        // x * Pr(H=1) = i - x * Pr(H=0)
        // x = i / (Pr(H=1) + Pr(H=0))

        let center = (count as f64 / (self.prob_call_exact + self.prob_call_mismatch)).round() as i32;
        (cmp::max(center - 30, 0) as u32, center as u32 + 30)
    }

    pub fn likelihood(&self, x: u32, count: u32, count_exact: u32) -> LogProb {
        let x = x;
        let count = count;
        let count_exact = count_exact;
        assert!(count >= count_exact);

        let mut summands = Vec::new();
        let probs = [self.prob_call_exact, self.prob_call_mismatch, self.prob_missed];
        for i in 0..(cmp::min(x, count) + 1) {
            let jmax = cmp::min(count_exact, i);
            let jmin = if count_exact + i > count { count_exact + i - count } else { 0 };
            for j in jmin..(jmax + 1) {
                let k = x - i;
                let mut p = multinomial_pdf(&probs, &[j, i - j, k]).ln() +
                        self.prob_miscall_exact *  (count_exact - j) as f64;
                let l = count - count_exact - (i - j);
                if l > 0 {
                    p += self.prob_miscall_mismatch * l as f64;
                }
                summands.push(p);
            }
        }
        let likelihood = logprobs::log_prob_sum(&summands);
        assert!(!likelihood.is_nan());
        likelihood
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use super::{MHD4, Params, Model};
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
        assert!(p.approx_eq(&0.21822058901379174));
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
}
