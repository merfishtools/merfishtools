use std::f64;
use std::collections;

use num::rational;

use bio::stats::logprobs::LogProb;
use bio::stats::logprobs;

use model;


pub type MeanExpression = rational::Ratio<u32>;
pub type PMF = collections::HashMap<MeanExpression, LogProb>;


pub fn pmf(expression_pmfs: &[model::expression::PMF]) -> PMF {
    let mut pmf = collections::HashMap::new();

    fn dfs(
        i: usize, expression_sum: u32, posterior_prob: LogProb,
        pmf: &mut collections::HashMap<MeanExpression, LogProb>,
        expression_pmfs: &[model::expression::PMF],
        min_prob: LogProb
    ) {
        if posterior_prob < min_prob {
            // stop if probability becomes too small
            return;
        }
        else if i < expression_pmfs.len() {
            let ref expr_pmf = expression_pmfs[i];
            let map = expr_pmf.map();
            for &(x, prob) in expr_pmf.iter() {
                dfs(
                    i + 1,
                    expression_sum + x,
                    posterior_prob + prob,
                    pmf,
                    expression_pmfs,
                    if x == map { f64::NEG_INFINITY } else { model::MIN_PROB }
                );
            }
        }
        else {
            let mean = rational::Ratio::new(expression_sum, expression_pmfs.len() as u32);
            if pmf.contains_key(&mean) {
                let p = pmf.get_mut(&mean).unwrap();
                *p = logprobs::log_prob_add(*p, posterior_prob);
            }
            else {
                pmf.insert(mean, posterior_prob);
            }
        }
    }

    dfs(0, 0, 0.0, &mut pmf, expression_pmfs, f64::NEG_INFINITY);

    pmf
}
