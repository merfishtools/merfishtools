#![allow(non_snake_case)]

use std::cmp;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs::{Prob, LogProb};


/// Probability to make exactly i 1->0 and j 0->1 errors.
fn xi(i: u8, j: u8, N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    p1.powi(i as i32) * p0.powi(j as i32) * (1.0 - p1).powi((m - i) as i32) * (1.0 - p0).powi((N - m - j) as i32)
}

/// Number of possibilities for i 1->0 and j 0->1 errors.
fn psi(i: u8, j: u8, N: u8, m: u8) -> f64 {
    combinations(m as u64, i as u64) * combinations((N - m)  as u64, j as u64)
}


/// Probability to see an exact readout given that we have a call.
fn prob_call_exact(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    xi(0, 0, N, m, p0, p1)
}


/// Probability to see a readout with one mismatch given that we have a call.
fn prob_call_mismatch(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    (m as f64 * xi(1, 0, N, m, p0, p1) + (N - m) as f64 * xi(0, 1, N, m, p0, p1)) / N as f64
}


/// Probability to see an exact readout given that we have a miscall.
fn prob_miscall_exact(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    xi(2, 2, N, m, p0, p1)
}


/// Probability to see a readout with one mismatch given that we have a miscall.
fn prob_miscall_mismatch(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    (psi(2, 1, N, m) * xi(2, 1, N, m, p0, p1) + psi(1, 2, N, m) * xi(1, 2, N, m, p0, p1) +
    psi(2, 3, N, m) * xi(2, 3, N, m, p0, p1) + psi(3, 2, N, m) * xi(3, 2, N, m, p0, p1)) /
    (psi(2, 1, N, m) + psi(1, 2, N, m) + psi(2, 3, N, m) + psi(3, 2, N, m))
}


/// Probability to completely miss a readout because of too many errors.
fn prob_missed(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    let mut denom = 0.0;
    let mut nom = 0.0;
    for k in 3..6 {
        for i in 0..cmp::min(m, k) + 1 {
            denom += psi(i, k - i, N, m) * xi(i, k - i, N, m, p0, p1);
            nom += psi(i, k - i, N, m);
        }
    }
    denom / nom
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
        Readout {
            prob_call_exact: prob_call_exact(N, m, p0, p1).ln(),
            prob_call_mismatch: prob_call_mismatch(N, m, p0, p1).ln(),
            prob_miscall_exact: prob_miscall_exact(N, m, p0, p1).ln(),
            prob_miscall_mismatch: prob_miscall_mismatch(N, m, p0, p1).ln(),
            prob_missed: prob_missed(N, m, p0, p1).ln()
        }
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use super::{prob_call_exact, prob_call_mismatch, prob_miscall_exact, prob_miscall_mismatch, prob_missed, xi};
    use bio::stats::logprobs::Prob;
    use nalgebra::ApproxEq;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;

    #[test]
    fn test_prob_call_exact() {
        let p = prob_call_exact(N, m, p0, p1);
        println!("{}", p);
        assert!(p.approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_prob_call_mismatch() {
        let p = prob_call_mismatch(N, m, p0, p1);
        println!("{}", p);
        assert!(p.approx_eq(&0.0237291000705869));
    }

    #[test]
    fn test_prob_miscall_exact() {
        let p = prob_miscall_exact(N, m, p0, p1);
        println!("{}", p);
        assert!(p.approx_eq(&0.000008616230962449854));
    }

    #[test]
    fn test_prob_miscall_mismatch() {
        let p = prob_miscall_mismatch(N, m, p0, p1);
        println!("{}", p);
        assert!(p.approx_eq(&0.00001879564966027473));
    }

    #[test]
    fn test_prob_missed() {
        let p = prob_missed(N, m, p0, p1);
        println!("{}", p);
        assert!(p.approx_eq(&0.000007861209464082392));
    }

    #[test]
    fn test_xi() {
        let p = xi(0, 0, N, m, p0, p1);
        assert!(p.approx_eq(&0.4019988717840602));
        let p = xi(1, 0, N, m, p0, p1);
        assert!(p.approx_eq(&0.04466654130934002));
        let p = xi(0, 1, N, m, p0, p1);
        assert!(p.approx_eq(&0.01674995299100251));
        let p = xi(2, 2, N, m, p0, p1);
        assert!(p.approx_eq(&8.616230962449852e-06));
        let p = xi(2, 1, N, m, p0, p1);
        assert!(p.approx_eq(&0.00020678954309879646));
        let p = xi(1, 2, N, m, p0, p1);
        assert!(p.approx_eq(&7.754607866204867e-05));
    }
}
