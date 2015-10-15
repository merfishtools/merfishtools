use std::cmp;
use std::f64;
use std::collections;

use itertools::Itertools;
use num::rational;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


pub type MeanExpression = rational::Ratio<u32>;


fn likelihood(x: u32, count: u32, count_exact: u32, readout_model: &model::Readout) -> LogProb {
    let x = x as u64;
    let count = count as u64;
    let count_exact = count_exact as u64;

    let x_c = cmp::min(count, x);
    let x_m = count - x_c;

    let imin = if count_exact + x_c > count { count_exact + x_c - count } else { 0 };
    let imax = cmp::min(count_exact, x) + 1;

    let summands = (imin..imax).map(|i| {
        let combs = combinations(count_exact, i).ln() +
                    combinations(count - count_exact, x_c - i).ln();
        let prob = readout_model.prob_call_exact * i as f64 +
                   readout_model.prob_call_mismatch * (x_c - i) as f64 +
                   readout_model.prob_miscall_exact * (count_exact - i) as f64 +
                   readout_model.prob_miscall_mismatch * (x_m + i - count_exact) as f64;
        (combs + prob).min(0.0)
    }).collect_vec();
    // removed (combinations(count, x_c).ln() as f64) since we already sum over all possibilities
    let likelihood = readout_model.prob_missed * (x - x_c) as f64 + logprobs::log_prob_sum(&summands).min(0.0);
    assert!(!likelihood.is_nan());
    likelihood
}


pub struct Expression {
    offset: u32,
    likelihoods: Vec<LogProb>,
    marginal: LogProb
}


impl Expression {
    pub fn new(count: u32, count_exact: u32, readout_model: &model::Readout) -> Self {
        let offset = if count_exact > 50 { count_exact - 50 } else { 0 };
        let likelihoods = (offset..count + 50).map(|x| likelihood(x, count, count_exact, readout_model)).collect_vec();
        // calculate (marginal / flat_prior)
        let marginal = logprobs::log_prob_sum(&likelihoods);
        // TODO trim values such that zero probabilities are not stored
        let mut expression = Expression {
            offset: offset,
            likelihoods: likelihoods,
            marginal: marginal
        };
        expression.refine_interval();
        expression
    }

    fn refine_interval(&mut self) {
        let min_prob = 0.000001f64.ln();
        let map = self.map();
        let mut min_x = self.min_x();
        for x in (self.min_x()..map).rev() {
            if self.posterior_prob(x) < min_prob {
                min_x = x;
                break;
            }
        }
        let mut max_x = self.max_x();
        for x in map..self.max_x() {
            if self.posterior_prob(x) < min_prob {
                max_x = x;
                break;
            }
        }
        self.likelihoods = self.likelihoods[(min_x - self.offset) as usize..(max_x - self.offset) as usize].to_vec();
        self.offset = min_x;
    }

    pub fn pmf(&self) -> Vec<(u32, LogProb)> {
        (self.min_x()..self.max_x()).map(|x| (x, self.posterior_prob(x))).collect_vec()
    }

    pub fn posterior_prob(&self, x: u32) -> LogProb {
        let x = x as usize;
        if x < self.offset as usize || x >= self.offset as usize + self.likelihoods.len() {
            f64::NEG_INFINITY
        }
        else {
            self.likelihoods[x - self.offset as usize] - self.marginal
        }
    }

    /// The maximum a posteriori probability estimate.
    pub fn map(&self) -> u32 {
        let mut max_x = self.min_x();
        let mut max = self.posterior_prob(max_x);
        for x in self.min_x() + 1..self.max_x() {
            let p = self.posterior_prob(x);
            if p >= max {
                max_x = x;
                max = p;
            }
        }
        max_x
    }

    pub fn expected_value(&self) -> f64 {
        (self.min_x()..self.max_x()).map(|x| x as f64 * self.posterior_prob(x).exp()).fold(0.0, |s, e| s + e)
    }

    pub fn variance(&self) -> f64 {
        (self.min_x()..self.max_x()).map(|x| (x as f64 - self.expected_value()).powi(2) * self.posterior_prob(x).exp()).fold(0.0, |s, e| s + e)
    }

    pub fn min_x(&self) -> u32 {
        self.offset
    }

    pub fn max_x(&self) -> u32 {
        self.offset + self.likelihoods.len() as u32
    }
}


pub struct ExpressionSet {
    pmf: collections::HashMap<MeanExpression, LogProb>
}


impl ExpressionSet {
    pub fn new(expressions: &[Expression]) -> Self {
        let mut pmf = collections::HashMap::new();

        fn dfs(
            i: usize, expression_sum: u32, posterior_prob: LogProb,
            pmf: &mut collections::HashMap<MeanExpression, LogProb>,
            expressions: &[Expression]
        ) {
            if i < expressions.len() {
                let ref expr = expressions[i];
                for x in expr.min_x()..expr.max_x() {
                    dfs(i + 1, expression_sum + x, posterior_prob + expr.posterior_prob(x), pmf, expressions);
                }
            }
            else {
                let mean = rational::Ratio::new(expression_sum, expressions.len() as u32);
                if pmf.contains_key(&mean) {
                    let p = pmf.get_mut(&mean).unwrap();
                    *p = logprobs::log_prob_add(*p, posterior_prob);
                }
                else {
                    pmf.insert(mean, posterior_prob);
                }
            }
        }

        dfs(0, 0, 0.0, &mut pmf, expressions);

        ExpressionSet {
            pmf: pmf
        }
    }

    /// Unorderered iteration over probability mass function (PMF).
    pub fn pmf(&self) -> collections::hash_map::Iter<MeanExpression, LogProb> {
        self.pmf.iter()
    }

    /// Conditional expectation of mean expression.
    pub fn expected_value(&self) -> f64 {
        self.pmf.iter().map(|(mean, prob)| {
            *mean.numer() as f64 / *mean.denom() as f64 * prob.exp()
        }).fold(0.0, |s, e| s + e)
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use itertools::Itertools;
    use nalgebra::ApproxEq;
    use bio::stats::logprobs::Prob;

    use super::{likelihood, Expression, ExpressionSet};
    use model;


    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;

    fn setup() -> model::Readout {
        model::Readout::new(N, m, p0, p1)
    }

    #[test]
    fn test_likelihood() {
        let readout = setup();
        assert!(likelihood(0, 0, 0, &readout).exp().approx_eq(&1.0));
        assert!(likelihood(1, 1, 1, &readout).exp().approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_posterior_prob() {
        let readout = setup();
        let expression = Expression::new(50, 10, &readout);
        println!("{:?}", (0..55).map(|x| expression.posterior_prob(x).exp()).collect_vec());
        assert!(false);
        let expression = Expression::new(5, 5, &readout);
        // check if x=5 yields highest probability
        assert_eq!((0..21).sorted_by(|&x, &y| {
            expression.posterior_prob(x).partial_cmp(&expression.posterior_prob(y)).unwrap()
        })[20], 5);
        // check that expression beyond window yields zero probability
        assert!(expression.posterior_prob(100).exp().approx_eq(&0.0));
    }

    #[test]
    fn test_set_expected_value() {
        let readout = setup();
        let expressions = [
            Expression::new(5, 1, &readout),
            Expression::new(5, 1, &readout),
            Expression::new(5, 1, &readout),
            Expression::new(1, 0, &readout)
        ];
        let expression_set = ExpressionSet::new(&expressions);
        // calculate expected value directly
        let expected_value = expressions.iter().map(|e| e.expected_value()).fold(0.0, |s, e| s + e) / expressions.len() as f64;
        assert!(expression_set.expected_value().approx_eq(&expected_value));
        // check if x=5 yields highest probability
        //assert_eq!((0..20).sorted_by(|&x, &y| {
        //    expression_set.posterior_prob(x).partial_cmp(&expression_set.posterior_prob(y)).unwrap()
        //})[19], 5);
    }
}
