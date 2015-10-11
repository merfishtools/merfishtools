use std::collections;
use std::cmp::Ordering;

use num::rational;
use itertools::Itertools;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


pub type LogFC = f64;
pub type FC = rational::Ratio<u32>;


pub struct Foldchange {
    pmf: collections::HashMap<FC, LogProb>
}


impl Foldchange {
    pub fn new(a: &model::ExpressionSet, b: &model::ExpressionSet) -> Self {
        let mut pmf = collections::HashMap::new();
        for ((a_mean, a_prob), (b_mean, b_prob)) in a.pmf().cartesian_product(b.pmf()) {
            let fc = b_mean / a_mean;
            let posterior_prob = b_prob + a_prob;

            if pmf.contains_key(&fc) {
                let p = pmf.get_mut(&fc).unwrap();
                *p = logprobs::log_prob_add(*p, posterior_prob);
            }
            else {
                pmf.insert(fc, posterior_prob);
            }
        }

        Foldchange {
            pmf: pmf
        }
    }

    /// Unorderered iteration over probability mass function (PMF).
    pub fn pmf(&self) -> collections::hash_map::Iter<FC, LogProb> {
        self.pmf.iter()
    }

    /// Minimum credible fold change f, i.e. Pr(F >= f | D) >= 0.95.
    // pub fn min_credible(&self) -> f64 {
    //     let sorted = self.pmf.iter().rev().sorted_by(|(a, _), (b, _)| a.partial_cmp(b));
    //     let cum_probs = logprobs::log_prob_cumsum(&sorted.iter().map(|(fc, prob)| prob).collect_vec());
    //
    // }

    /// Conditional expectation of fold change.
    pub fn expected_value(&self) -> f64 {
        self.pmf.iter().map(|(fc, prob)| {
            *fc.numer() as f64 / *fc.denom() as f64 * prob.exp()
        }).sum()
    }

    /// Conditional variance of fold change.
    pub fn variance(&self) -> f64 {
        let expected_value = self.expected_value();
        self.pmf.iter().map(|(fc, prob)| {
            (*fc.numer() as f64 / *fc.denom() as f64 - expected_value) * prob.exp()
        }).sum()
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use bio::stats::logprobs::Prob;
    use nalgebra::ApproxEq;

    use model;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;



    fn setup() -> model::Readout {
        model::Readout::new(N, m, p0, p1)
    }

    #[test]
    fn test_expected_value() {
        let readout = setup();
        let a = [
            model::Expression::new(5, 5, &readout),
            model::Expression::new(5, 5, &readout),
            model::Expression::new(5, 5, &readout),
            model::Expression::new(5, 5, &readout)
        ];

        let b = [
            model::Expression::new(15, 15, &readout),
            model::Expression::new(15, 15, &readout),
            model::Expression::new(15, 15, &readout),
            model::Expression::new(15, 15, &readout)
        ];

        let a = model::ExpressionSet::new(&a);
        let b = model::ExpressionSet::new(&b);

        let fc = Foldchange::new(&a, &b);

        println!("{}", fc.expected_value());
        assert!(fc.expected_value().approx_eq(&2.999373414428798));
    }
}
