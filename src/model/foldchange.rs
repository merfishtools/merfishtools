use std::collections;
use std::iter;

use num::rational;
use itertools::Itertools;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;
use model::pmf::PMF;


pub type LogFC = f64;
pub type FC = rational::Ratio<u32>;


pub struct Foldchange {
    inner: collections::HashMap<FC, LogProb>
}


impl<'a> PMF<'a, FC, iter::Cloned<collections::hash_map::Keys<'a, FC, LogProb>>> for Foldchange {
    fn domain(&'a self) -> iter::Cloned<collections::hash_map::Keys<'a, FC, LogProb>> {
        self.inner.keys().cloned()
    }

    fn posterior_prob(&'a self, fc: FC) -> LogProb {
        *self.inner.get(&fc).unwrap()
    }

    fn cast(value: FC) -> f64 {
        *value.numer() as f64 / *value.denom() as f64
    }
}


impl Foldchange {
    pub fn new(a: &model::ExpressionSet, b: &model::ExpressionSet) -> Self {
        let mut pmf = collections::HashMap::new();
        for (a_mean, b_mean) in a.domain().cartesian_product(b.domain()) {
            let fc = b_mean / a_mean;
            let posterior_prob = b.posterior_prob(b_mean) + a.posterior_prob(a_mean);

            if pmf.contains_key(&fc) {
                let p = pmf.get_mut(&fc).unwrap();
                *p = logprobs::log_prob_add(*p, posterior_prob);
            }
            else {
                pmf.insert(fc, posterior_prob);
            }
        }

        Foldchange {
            inner: pmf
        }
    }
/*
    /// Probability mass function (PMF).
    pub fn pmf(&self) -> PMF {
        let mut pmf = self.pmf.iter().map(|(fc, prob)| (ratio_as_f64(fc).log2(), *prob)).collect_vec();
        pmf.sort_by(|&(a, _), &(b, _)| a.partial_cmp(&b).unwrap());
        model::pmf::PMF::new(pmf)
    }*/
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

    // #[test]
    // fn test_expected_value() {
    //     let readout = setup();
    //     let a = [
    //         model::Expression::new(5, 5, &readout),
    //         model::Expression::new(5, 5, &readout),
    //         model::Expression::new(5, 5, &readout),
    //         model::Expression::new(5, 5, &readout)
    //     ];
    //
    //     let b = [
    //         model::Expression::new(15, 15, &readout),
    //         model::Expression::new(15, 15, &readout),
    //         model::Expression::new(15, 15, &readout),
    //         model::Expression::new(15, 15, &readout)
    //     ];
    //
    //     let a = model::ExpressionSet::new(&a);
    //     let b = model::ExpressionSet::new(&b);
    //
    //     let fc = Foldchange::new(&a, &b);
    //
    //     println!("{}", fc.expected_value());
    //     assert!(fc.expected_value().approx_eq(&2.999373414428798));
    // }
}
