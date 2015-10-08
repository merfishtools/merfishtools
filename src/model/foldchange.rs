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
                self.a.posterior_prob(x) + self.b.posterior_prob((f * x as f64) as u32)
            }
            else {
                f64::NEG_INFINITY
            }
        }).collect_vec();
        logprobs::log_prob_sum(&p)
    }
/*
    pub fn expectation(&self) -> f64 {
        let b_prob = (self.b.min_x()..self.b.max_x() + 1).map(|x| self.b.posterior_prob(x)).collect_vec();

        for x_a in (self.a.min_x()..self.a.max_x() + 1) {
            let a_prob = self.a.posterior_prob(x_a);
            for (i, x_b) in (self.b.min_x()..self.b.max_x() + 1).enumerate() {

            }
        }
    }*/
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use bio::stats::logprobs::Prob;

    use model;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;



    fn setup() -> model::Readout {
        model::Readout::new(N, m, p0, p1)
    }

    /*#[test]
    fn test_posterior_prob() {
        let readout = setup();
        let mut a = model::ExpressionSet::new();
        a.push(model::Expression::new(5, 5, &readout));
        a.push(model::Expression::new(5, 5, &readout));
        a.push(model::Expression::new(5, 5, &readout));
        a.push(model::Expression::new(5, 5, &readout));

        let mut b = model::ExpressionSet::new();
        b.push(model::Expression::new(15, 15, &readout));
        b.push(model::Expression::new(15, 15, &readout));
        b.push(model::Expression::new(15, 15, &readout));
        b.push(model::Expression::new(15, 15, &readout));

        let fc = Foldchange::new(&a, &b);

        println!("{}", fc.posterior_prob(3.0));
        assert!(false);
    }*/
}
