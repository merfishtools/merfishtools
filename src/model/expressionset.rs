use std::collections;
use std::collections::hash_map::Entry::{Occupied, Vacant};
use std::mem;

use num::rational;
use itertools::Itertools;

use bio::stats::logprobs;

use model;


pub type MeanExpression = rational::Ratio<u32>;
pub type PMF = model::pmf::PMF<MeanExpression>;


fn norm_scale_factors(scale_factors: &[f64]) -> Vec<rational::Ratio<u32>> {
    let mut factors = Vec::new();
    let min = scale_factors.iter().fold(1.0, |m, &e| if e < m { e } else { m });
    for &s in scale_factors.iter() {
        let s = ((s / min) * 100.0).round() / 100.0;
        factors.push(rational::Ratio::new(s.trunc() as u32 * 100 + s.fract() as u32, 100))
    }

    factors
}


pub fn pmf(expression_pmfs: &[model::expression::PMF], scale_factors: &[f64]) -> PMF {
    let scale_factors = norm_scale_factors(scale_factors);
    println!("{:?}", scale_factors);

    let mut curr = collections::HashMap::new();
    let mut prev = collections::HashMap::new();
    curr.insert(rational::Ratio::from_integer(0), 0.0);

    for (pmf, scale) in expression_pmfs.iter().zip(scale_factors) {
        mem::swap(&mut curr, &mut prev); // TODO does this really swap curr and prev?

        curr.clear();

        for (s, p) in prev.iter() {
            for &(x, x_prob) in pmf.iter() {
                let s = s + rational::Ratio::from_integer(x) * scale;
                match curr.entry(s) {
                    Occupied(mut entry) => {
                        let prob = logprobs::log_prob_add(*entry.get(), p + x_prob);
                        entry.insert(prob);
                    },
                    Vacant(entry)   => {
                        entry.insert(p + x_prob);
                    }
                }
            }
        }
    }

    PMF::new(curr.iter().filter_map(|(s, p)| {
        if *p >= model::MIN_PROB {
            Some((s / rational::Ratio::from_integer(expression_pmfs.len() as u32), *p))
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
        let pmf = pmf(&pmfs, &[1.0, 1.0, 1.0, 1.0]);

        let total = log_prob_sum(&pmf.iter().map(|&(_, prob)| prob).collect_vec());
        let values = pmf.iter().map(|e| (*e.0.numer() as f64 / *e.0.denom() as f64, e.1)).collect_vec();

        println!("{:?}", values);
        println!("{:?}", total);
        assert!(false);
        assert!(total.approx_eq(&-0.000003372325827477596));
    }
}
