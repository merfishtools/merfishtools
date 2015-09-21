use std::f64;

use itertools::Itertools;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


pub struct Foldchange<'a> {
    a: &'a model::ExpressionSet,
    b: &'a model::ExpressionSet
}


impl<'a> Foldchange<'a> {
    pub fn new(a: &'a model::ExpressionSet, b: &'a model::ExpressionSet) -> Self {
        Foldchange {
            a: a,
            b: b
        }
    }

    pub fn posterior_prob(&self, f: f64) -> LogProb {
        let xmin = self.a.min_x();
        let xmax = self.a.max_x();

        let p = (xmin..xmax + 1).map(|x| {
            let x_ = x as f64 * f;
            if x_ % 1.0 == 0.0 {
                // multiply in log space
                self.a.posterior_prob(x) + self.b.posterior_prob((f * x as f64) as usize)
            }
            else {
                f64::NEG_INFINITY
            }
        }).collect_vec();
        logprobs::log_prob_sum(&p)
    }
}
