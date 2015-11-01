use std::f64;

use num::rational;
use itertools::Itertools;

use bio::stats::logprobs::LogProb;
use bio::stats::logprobs;

use model;


pub type MeanExpression = rational::Ratio<u32>;
pub type PMF = model::pmf::PMF<MeanExpression>;


pub fn pmf(expression_pmfs: &[model::expression::PMF]) -> PMF {
    let max_sum = expression_pmfs.iter()
                                 .map(|pmf| pmf.iter().last().unwrap().0)
                                 .fold(0, |s, e| s + e) as usize;

    let mut probs = [vec![f64::NEG_INFINITY; max_sum + 1], vec![f64::NEG_INFINITY; max_sum + 1]];
    probs[1][0] = 0.0;
    let mut curr = 0;

    for (i, pmf) in expression_pmfs.iter().enumerate() {
        curr = i % 2;
        let prev = 1 - curr;

        for p in probs[curr].iter_mut() {
            *p = f64::NEG_INFINITY;
        }

        for s in 0..max_sum {
            let p = probs[prev][s];
            if p != f64::NEG_INFINITY {
                for &(x, x_prob) in pmf.iter() {
                    let s = s + x as usize;
                    if s <= max_sum {
                        let prob = &mut probs[curr][s];
                        *prob = logprobs::log_prob_add(*prob, p + x_prob);
                    }
                }
            }
        }
    }

    PMF::new(probs[curr].iter().enumerate().filter_map(|(s, p)| {
        if *p >= model::MIN_PROB {
            Some((rational::Ratio::new(s as u32, expression_pmfs.len() as u32), *p))
        }
        else {
            None
        }
    }).collect_vec())
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use itertools::Itertools;
    use nalgebra::ApproxEq;
    use bio::stats::logprobs::{Prob, log_prob_sum};


    use super::*;
    use model;


    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;

    fn setup() -> model::Readout {
        model::Readout::new(N, m, p0, p1)
    }

    #[test]
    fn test_pmf() {
        let readout = setup();
        let pmfs = [
            model::expression::pmf(5, 1, &readout),
            model::expression::pmf(10, 1, &readout),
            model::expression::pmf(3, 1, &readout),
            model::expression::pmf(24, 1, &readout)
        ];
        let pmf = pmf(&pmfs);

        let total = log_prob_sum(&pmf.iter().map(|&(_, prob)| prob).collect_vec());

        println!("{:?}", total);
        assert!(total.approx_eq(&0.0));
    }
}
