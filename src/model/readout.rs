#![allow(non_snake_case)]

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs::{Prob, LogProb};


/// Probability to see an exact readout given that we have a call.
fn prob_call_exact(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    (1.0 - p1).powi(m as i32) * (1.0 - p0).powi((N - m) as i32)
}


/// Probability to see a readout with one mismatch given that we have a call.
fn prob_call_mismatch(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    m as Prob * p1 * (1.0 - p1).powi((m - 1) as i32) * (1.0 - p0).powi((N - m) as i32) +
    (N - m) as Prob * p0 * (1.0 - p1).powi(m as i32) * (1.0 - p0).powi((N - m - 1) as i32)
}


/// Helper function.
fn xi(k: u8, N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    let mut xi = 0.0;
    for i in 0..k {
        xi += combinations(m as u64, i as u64) * combinations((N - m) as u64, (k - i) as u64) *
              p1.powi(i as i32) * p0.powi((k - i) as i32) *
              (1.0 - p1).powi((m - i) as i32) * (1.0 - p0).powi((N - m - k + i) as i32);
    }
    xi
}


/// Probability to see an exact readout given that we have a call.
fn prob_miscall_exact(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    xi(4, N, m, p0, p1)
}


/// Probability to see a readout with one mismatch given that we have a miscall.
fn prob_miscall_mismatch(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    xi(3, N, m, p0, p1) + xi(5, N, m, p0, p1)
}


/// Probability to miss a readout because it is assigned to another gene.
fn prob_missed(N: u8, m: u8, p0: Prob, p1: Prob) -> Prob {
    xi(3, N, m, p0, p1) + xi(4, N, m, p0, p1) + xi(5, N, m, p0, p1)
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

    use super::{prob_call_exact, prob_call_mismatch, xi};
    use bio::stats::logprobs::Prob;
    use nalgebra::ApproxEq;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;

    #[test]
    fn test_prob_call_exact() {
        let p = prob_call_exact(N, m, p0, p1);
        assert!(p.approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_prob_call_mismatch() {
        let p = prob_call_mismatch(N, m, p0, p1);
        assert!(p.approx_eq(&0.3796656011293902));
    }

    #[test]
    fn test_xi() {
        assert!(xi(3, N, m, p0, p1).approx_eq(&0.041758563359513209));
        assert!(xi(4, N, m, p0, p1).approx_eq(&0.0079580316940839561));
        assert!(xi(5, N, m, p0, p1).approx_eq(&0.0010638203079106705));
    }
}
