use std::cmp;
use std::f64;

use rand;
use rand::distributions::IndependentSample;
use itertools::Itertools;
use itertools;

use bio::stats::combinatorics::combinations;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


fn likelihood(x: u64, count: u64, count_exact: u64, readout_model: &model::Readout) -> LogProb {
    let x_c = cmp::min(count, x);
    let x_m = count - x_c;

    let imin = if count_exact + x_c > count { count_exact + x_c - count } else { 0 };
    let imax = cmp::min(count_exact, x) + 1;
    // work in log-space
    let summands = (imin..imax).map(|i| {
        combinations(count_exact, i).ln() +
        combinations(count - count_exact, x_c - i).ln() +
        readout_model.prob_call_exact * i as f64 +
        readout_model.prob_call_mismatch * (x_c - i) as f64 +
        readout_model.prob_miscall_exact * (count_exact - i) as f64 +
        readout_model.prob_miscall_mismatch * (x_m + i - count_exact) as f64
    }).collect_vec();
    // TODO think about the combinations (divide or multiply??)
    let likelihood = readout_model.prob_missed * (x - x_c) as f64 + combinations(count, x_c).ln() as f64 +
                     logprobs::log_prob_sum(&summands);
    likelihood
}


pub struct Expression {
    offset: usize,
    likelihoods: Vec<LogProb>,
    marginal: LogProb
}


impl Expression {
    pub fn new(count: usize, count_exact: usize, readout_model: &model::Readout) -> Self {
        let offset = if count_exact > 50 { count_exact - 50 } else { 0 };
        let likelihoods = (offset..count + 50).map(|x| likelihood(x as u64, count as u64, count_exact as u64, readout_model)).collect_vec();
        // calculate (marginal / flat_prior)
        let marginal = logprobs::log_prob_sum(&likelihoods);
        Expression {
            offset: offset,
            likelihoods: likelihoods,
            marginal: marginal
        }
    }

    pub fn posterior_prob(&self, x: usize) -> LogProb {
        if x < self.offset || x >= self.offset + self.likelihoods.len() {
            f64::NEG_INFINITY
        }
        else {
            self.likelihoods[x - self.offset] - self.marginal
        }
    }

    pub fn min_x(&self) -> usize {
        self.offset
    }

    pub fn max_x(&self) -> usize {
        self.offset +  self.likelihoods.len()
    }
}


pub struct ExpressionSet {
    expressions: Vec<Expression>
}


impl ExpressionSet {
    pub fn new() -> Self {
        ExpressionSet {
            expressions: Vec::new()
        }
    }

    pub fn push(&mut self, expression: Expression) {
        self.expressions.push(expression);
    }

    pub fn posterior_prob(&self, x_expectation: usize) -> LogProb {
        if self.expressions.len() == 1 {
            self.expressions[0].posterior_prob(x_expectation)
        }
        else {
            let posteriors = self.expressions.iter().map(|e| e.posterior_prob(x_expectation)).collect_vec();
            let n = posteriors.len();

            let mut rng = rand::thread_rng();
            let idx = rand::distributions::Range::new(0, n);
            // bootstrap 500 times
            let bootstrap_n = 500;
            let bootstrap_values = itertools::RepeatCall::new(|| {
                itertools::RepeatCall::new(|| posteriors[idx.ind_sample(&mut rng)]).take(n).sum()
            }).take(bootstrap_n).collect_vec();
            // calculate the expected posterior probability over the bootstap
            logprobs::log_prob_sum(&bootstrap_values) - (bootstrap_n as f64).ln()
        }
    }

    pub fn min_x(&self) -> usize {
        // TODO use something more sophisticated here (e.g. credible intervals)
        self.expressions.iter().map(|e| e.min_x()).min().unwrap()
    }

    pub fn max_x(&self) -> usize {
        // TODO use something more sophisticated here (e.g. credible intervals)
        self.expressions.iter().map(|e| e.max_x()).max().unwrap()
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use test::Bencher;

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
        let expression = Expression::new(50, 30, &readout);
        println!("{:?}", (0..51).map(|x| expression.posterior_prob(x).exp()).collect_vec());
        //assert!(false);
        let expression = Expression::new(5, 5, &readout);
        // check if x=5 yields highest probability
        assert_eq!((0..21).sorted_by(|&x, &y| {
            expression.posterior_prob(x).partial_cmp(&expression.posterior_prob(y)).unwrap()
        })[20], 5);
        // check that expression beyond window yields zero probability
        assert!(expression.posterior_prob(100).exp().approx_eq(&0.0));
    }

    #[test]
    fn test_set_posterior_prob() {
        let readout = setup();
        let mut expression_set = ExpressionSet::new();
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(1, 0, &readout));
        // check if x=5 yields highest probability
        assert_eq!((0..20).sorted_by(|&x, &y| {
            expression_set.posterior_prob(x).partial_cmp(&expression_set.posterior_prob(y)).unwrap()
        })[19], 5);
    }

    #[bench]
    fn bench_set_posterior_prob(b: &mut Bencher) {
        let readout = setup();
        let mut expression_set = ExpressionSet::new();
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(5, 5, &readout));
        expression_set.push(Expression::new(1, 0, &readout));

        b.iter(|| expression_set.posterior_prob(5));
    }
}
