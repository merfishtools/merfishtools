use std::mem;

use itertools::Itertools;
use ndarray::prelude::*;

use bio::stats::{LogProb, Prob};

use crate::io::codebook::Codeword;

/// Basic model for the probability of making i 1-0 and j 0-1 errors.
pub struct Xi {
    p0: Vec<LogProb>,
    p1: Vec<LogProb>,
}

impl Xi {
    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `p0` - error probabilities for 0-1 errors
    /// * `p1` - error probabilities for 1-0 errors
    pub fn new(p0: &[Prob], p1: &[Prob]) -> Self {
        let tolog = |&p| LogProb::from(p);
        Xi {
            p0: p0.iter().map(&tolog).collect_vec(),
            p1: p1.iter().map(&tolog).collect_vec(),
        }
    }

    /// Calculate the probability of transforming source into target during a MERFISH experiment.
    pub fn prob(&self, source: &Codeword, target: &Codeword) -> [LogProb; 2] {
        let mut curr = Array1::from_elem(3, LogProb::ln_zero());
        let mut prev = Array1::from_elem(3, LogProb::ln_zero());
        prev[1] = LogProb::ln_one();

        for k in 0..source.len() {
            for d in 1..curr.shape()[0] {
                let (d_0, d_1, p_0) = match (source.get(k).unwrap(), target.get(k).unwrap()) {
                    (true, true) => (d - 1, d, self.p1[k]),
                    (false, false) => (d - 1, d, self.p0[k]),
                    (false, true) => (d, d - 1, self.p0[k]),
                    (true, false) => (d, d - 1, self.p1[k]),
                };

                let p = (p_0 + prev[d_0]).ln_add_exp(p_0.ln_one_minus_exp() + prev[d_1]);

                curr[d] = p;
            }
            mem::swap(&mut prev, &mut curr);
        }
        [prev[1], prev[2]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio::stats::{LogProb, Prob};
    use bit_vec::BitVec;

    #[test]
    fn test_xi_exact_or_mismatch() {
        let a = BitVec::from_bytes(&[0b10101101]);
        let b = BitVec::from_bytes(&[0b10101011]);
        let p0 = LogProb::from(Prob(0.005));
        let p1 = LogProb::from(Prob(0.01));

        let xi = Xi::new(&[Prob::from(p0); 8], &[Prob::from(p1); 8]);

        let p = xi.prob(&a, &b);

        let truth =
            LogProb(*p1.ln_one_minus_exp() * 4.0) + LogProb(*p0.ln_one_minus_exp() * 2.0) + p1 + p0;

        assert_relative_eq!(*p[0], *truth);

        let truth_mismatch = LogProb::ln_sum_exp(&[
            // mismatch in matching 1-bits, 4 positions
            LogProb::from(Prob(*Prob::from(truth - p1.ln_one_minus_exp() + p1) * 4.0)),
            // mismatch in matching 0-bits, 2 positions
            LogProb::from(Prob(*Prob::from(truth - p0.ln_one_minus_exp() + p0) * 2.0)),
            // mismatch in mismatching 1-bits, 1 position
            truth - p1 + p1.ln_one_minus_exp(),
            // mismatch in mismatching 0-bits, 1 position
            truth - p0 + p0.ln_one_minus_exp(),
        ]);

        assert_relative_eq!(*p[1], *truth_mismatch, epsilon = 0.001);
    }

    #[test]
    fn test_xi_exact_total() {
        let p0 = LogProb::from(Prob(0.005));
        let p1 = LogProb::from(Prob(0.01));
        let xi = Xi::new(&[Prob::from(p0); 8], &[Prob::from(p1); 8]);

        let a = BitVec::from_bytes(&[0b10101101]);

        let mut probs = Vec::new();
        for b in 0..=255 {
            let b = BitVec::from_bytes(&[b]);
            probs.push(xi.prob(&a, &b)[0]);
        }
        let p = LogProb::ln_sum_exp(&probs);

        assert_relative_eq!(*p, *LogProb::ln_one());
    }

    #[test]
    fn test_xi_noise() {
        let p0 = LogProb::from(Prob(0.04));
        let p1 = LogProb::from(Prob(0.1));
        let xi = Xi::new(&[Prob::from(p0); 16], &[Prob::from(p1); 16]);

        let noise = BitVec::from_elem(16, false);
        let target = BitVec::from_bytes(&[0b10001100, 0b00000010]);

        let p = xi.prob(&noise, &target);

        let truth_exact = LogProb(*p0 * 4.0) + LogProb(*p0.ln_one_minus_exp() * 12.0);

        let truth_mismatch = LogProb::ln_sum_exp(&[
            // mismatch in matching 0-bits, 12 possible positions
            LogProb::from(Prob(
                *Prob::from(truth_exact - p0.ln_one_minus_exp() + p0) * 12.0,
            )),
            // mismatch in mismatching 1-bits
            LogProb::from(Prob(
                *Prob::from(truth_exact - p0 + p0.ln_one_minus_exp()) * 4.0,
            )),
        ]);

        assert_relative_eq!(*p[0], *truth_exact, epsilon = 0.0001);
        assert_relative_eq!(*p[1], *truth_mismatch);
    }
}
