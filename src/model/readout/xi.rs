use std::mem;

use itertools::Itertools;
use ndarray::prelude::*;
use bit_vec::BitVec;

use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codeword};


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

    /// Return matrix of error probabilities.
    /// Entry (i, j) contains the probability for i 1-0 and j 0-1 errors.
    ///
    /// # Arguments
    /// * `max_err` - maximum number of errors considered
    pub fn prob_error(
        &self,
        source: &codebook::Record,
        target: Option<&codebook::Record>,
        max_err: usize
    ) -> Array2<LogProb> {
        let mask = if let Some(target) = target {
            Some(source.diff(target))
        } else {
            None
        };

        self.calc(source.codeword(), mask.as_ref(), max_err)
    }

    fn calc(&self, codeword: &Codeword, mask: Option<&BitVec>, max_err: usize) -> Array2<LogProb> {
        let n = max_err + 1;
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
