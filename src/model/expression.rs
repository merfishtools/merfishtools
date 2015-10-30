use std::f64;
use std::ops::Range;

use itertools::Itertools;

use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model::pmf::PMF;
use model;


pub struct Expression {
    offset: u32,
    likelihoods: Vec<LogProb>,
    marginal: LogProb
}


impl<'a> PMF<'a, u32, Range<u32>> for Expression {
    fn domain(&'a self) -> Range<u32> {
        self.min_x()..self.max_x()
    }

    fn posterior_prob(&'a self, x: u32) -> LogProb {
        let x = x as usize;
        if x < self.offset as usize || x >= self.offset as usize + self.likelihoods.len() {
            f64::NEG_INFINITY
        }
        else {
            self.likelihoods[x - self.offset as usize] - self.marginal
        }
    }

    fn cast(value: u32) -> f64 {
        value as f64
    }
}


impl Expression {
    pub fn new(count: u32, count_exact: u32, readout_model: &model::Readout) -> Self {
        let offset = if count_exact > 50 { count_exact - 50 } else { 0 };
        let likelihoods = (offset..count + 50).map(|x| readout_model.likelihood(x, count, count_exact)).collect_vec();
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
        let map = self.map();
        let mut min_x = self.min_x();
        for x in (self.min_x()..map).rev() {
            if self.posterior_prob(x) < model::MIN_PROB {
                min_x = x;
                break;
            }
        }
        let mut max_x = self.max_x();
        for x in map..self.max_x() {
            if self.posterior_prob(x) < model::MIN_PROB {
                max_x = x;
                break;
            }
        }
        self.likelihoods = self.likelihoods[(min_x - self.offset) as usize..(max_x - self.offset) as usize].to_vec();
        self.offset = min_x;
    }

    /*
    pub fn pmf(&self) -> PMF {
        model::pmf::PMF::new((self.min_x()..self.max_x()).map(|x| (x, self.posterior_prob(x))).collect_vec())
    }*/

/*
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
    }*/

    pub fn min_x(&self) -> u32 {
        self.offset
    }

    pub fn max_x(&self) -> u32 {
        self.offset + self.likelihoods.len() as u32
    }
}

/*
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
        let expression = Expression::new(100, 99, &readout);
        println!("{:?}", (90..110).map(|x| expression.posterior_prob(x).exp()).collect_vec());
        assert!(false);
        let expression = Expression::new(5, 5, &readout);
        // check if x=5 yields highest probability
        assert_eq!((0..21).sorted_by(|&x, &y| {
            expression.posterior_prob(x).partial_cmp(&expression.posterior_prob(y)).unwrap()
        })[20], 5);
        // check that expression beyond window yields zero probability
        assert!(expression.posterior_prob(100).exp().approx_eq(&0.0));
    }

    // #[test]
    // fn test_set_expected_value() {
    //     let readout = setup();
    //     let expressions = [
    //         Expression::new(5, 1, &readout),
    //         Expression::new(5, 1, &readout),
    //         Expression::new(5, 1, &readout),
    //         Expression::new(1, 0, &readout)
    //     ];
    //     let expression_set = ExpressionSet::new(&expressions);
    //     // calculate expected value directly
    //     let expected_value = expressions.iter().map(|e| e.expected_value()).fold(0.0, |s, e| s + e) / expressions.len() as f64;
    //     assert!(expression_set.expected_value().approx_eq(&expected_value));
    // }
}
*/
