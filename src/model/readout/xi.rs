use std::mem;

use itertools::Itertools;
use ndarray::prelude::*;

use bio::stats::{Prob, LogProb};

use io::codebook::Codeword;


/// Basic model for the probability of making i 1-0 and j 0-1 errors.
pub struct Xi {
    p0: Vec<LogProb>,
    p1: Vec<LogProb>
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
            p1: p1.iter().map(&tolog).collect_vec()
        }
    }

    pub fn prob(&self, source: &Codeword, target: &Codeword) -> [LogProb; 2] {
        let mut curr = Array1::from_elem(3, LogProb::ln_zero());
        let mut prev = Array1::from_elem(3, LogProb::ln_zero());
        prev[1] = LogProb::ln_one();

        for k in 0..source.len() {
            for d in 1..curr.shape()[0] {
                let (d_0, d_1, p_0) = match (source.get(k).unwrap(), target.get(k).unwrap()) {
                    (true, true)   => (d - 1, d, self.p1[k]),
                    (false, false) => (d - 1, d, self.p0[k]),
                    (false, true)  => (d, d - 1, self.p0[k]),
                    (true, false)  => (d, d - 1, self.p1[k])
                };

                let p = (p_0 + prev[d_0]).ln_add_exp(p_0.ln_one_minus_exp() + prev[d_1]);

                curr[d] = p;
            }
            mem::swap(&mut prev, &mut curr);
        }
        [prev[1], prev[2]]
    }
}
