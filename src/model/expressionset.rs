use std::mem;
use std::f64;

use num::rational;
use itertools::Itertools;

use bio::stats::logprobs;

use model;


pub type MeanExpression = rational::Ratio<u32>;
pub type PMF = model::pmf::PMF<MeanExpression>;

const SCALE: f64 = 10.0;

pub fn pmf(expression_pmfs: &[model::expression::PMF]) -> PMF {
    let max_sum = expression_pmfs.iter()
                                 .map(|pmf| pmf.iter().last().unwrap().value)
                                 .fold(0, |s, e| s + (e * SCALE).round() as usize);

    let mut curr = vec![f64::NEG_INFINITY; max_sum + 1];
    let mut prev = vec![f64::NEG_INFINITY; max_sum + 1];
    curr[0] = 0.0;

    for pmf in expression_pmfs.iter() {
        mem::swap(&mut curr, &mut prev);

        for p in curr.iter_mut() {
            *p = f64::NEG_INFINITY;
        }

        for s in 0..max_sum {
            let p = prev[s];
            if p != f64::NEG_INFINITY {
                for x in pmf.iter() {
                    let s = s + (x.value * SCALE) as usize;
                    if s <= max_sum {
                        let prob = &mut curr[s];
                        *prob = logprobs::log_prob_add(*prob, p + x.prob);
                    }
                }
            }
        }
    }

    PMF::new(curr.iter().enumerate().filter_map(|(s, p)| {
        if *p >= model::MIN_PROB {
            Some(model::pmf::Entry{
                value: rational::Ratio::new(s as u32, SCALE as u32 * expression_pmfs.len() as u32),
                prob: *p
            })
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

        let total = log_prob_sum(&pmf.iter().map(|e| e.prob).collect_vec());
        let values = pmf.iter().map(|e| (*e.value.numer() as f64 / *e.value.denom() as f64, e.prob)).collect_vec();

        println!("{:?}", values);
        println!("{:?}", total);
        assert!(total.approx_eq(&-0.000003372325827477596));
    }
}
