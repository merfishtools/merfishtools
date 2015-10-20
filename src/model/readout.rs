#![allow(non_snake_case)]

use std::cmp;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs::{Prob, LogProb};


struct Factory {
    N: u8,
    m: u8,
    p0: Prob,
    p1: Prob
}


impl Factory {

    /// Probability to make exactly i 1->0 and j 0->1 errors.
    fn xi(&self, i: u8, j: u8) -> Prob {
        self.p1.powi(i as i32) * self.p0.powi(j as i32) * (1.0 - self.p1).powi((self.m - i) as i32) * (1.0 - self.p0).powi((self.N - self.m - j) as i32)
    }

    /// Number of possibilities for i 1->0 and j 0->1 errors.
    fn psi(&self, i: u8, j: u8) -> f64 {
        combinations(self.m as u64, i as u64) * combinations((self.N - self.m)  as u64, j as u64)
    }


    /// Probability to see an exact readout given that we have a call.
    fn prob_call_exact(&self) -> Prob {
        self.xi(0, 0)
    }


    /// Probability to see a readout with one mismatch given that we have a call.
    fn prob_call_mismatch(&self) -> Prob {
        //(self.m as f64 * self.xi(1, 0) + (self.N - self.m) as f64 * self.xi(0, 1)) / self.N as f64
        (self.m as f64 * self.xi(1, 0) + (self.N - self.m) as f64 * self.xi(0, 1))
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
        for k in 3..6 {
            for i in 0..cmp::min(self.m, k) + 1 {
                p += self.psi(i, k - i) * self.xi(i, k - i);
            }
        }
        p
    }
}


/// Readout probabilities.
pub struct Readout {
    pub prob_call_exact: LogProb,
    pub prob_call_mismatch: LogProb,
    pub prob_miscall_exact: LogProb,
    pub prob_miscall_mismatch: LogProb,
    pub prob_missed: LogProb
}


impl Readout {
    pub fn new(N: u8, m: u8, p0: Prob, p1: Prob) -> Self {
        let factory = Factory { N: N, m: m, p0: p0, p1: p1 };
        Readout {
            prob_call_exact: factory.prob_call_exact().ln(),
            prob_call_mismatch: factory.prob_call_mismatch().ln(),
            prob_miscall_exact: factory.prob_miscall_exact().ln(),
            prob_miscall_mismatch: factory.prob_miscall_mismatch().ln(),
            prob_missed: factory.prob_missed().ln()
        }
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use super::{Factory};
    use bio::stats::logprobs::Prob;
    use nalgebra::ApproxEq;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;
    const factory: Factory = Factory { N: 16, m: 4, p0: 0.04, p1: 0.1 };


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
        assert!(p.approx_eq(&0.053047441463627984));
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
