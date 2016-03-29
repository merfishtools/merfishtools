use std::collections;

use num::traits::ToPrimitive;
use num::rational;
use itertools::Itertools;
use bio::stats::logprobs;

use model;


pub type CV = f64;


#[allow(unused_parens)]
pub fn pmf(pmfs: &[model::expressionset::PMF]) -> model::diffexp::PMF {
    let mut pmf = collections::HashMap::new();
    let combinations = iproduct!(pmfs);
    let n = pmfs.len() as f64;

    for groups in combinations {
        let prob = groups.iter().map(|group| group.prob).fold(0.0, |s, p| s + p);
        if prob >= model::MIN_PROB {
            // sample mean
            let mean = groups.iter().map(|group| {
                *group.value.numer() as f64 / *group.value.denom() as f64
            }).fold(0.0, |s, e| s + e) / n;
            // sample standard deviation
            let std = (groups.iter().map(|group| {
                    (*group.value.numer() as f64 / *group.value.denom() as f64 - mean).powi(2)
            }).fold(0.0, |s, e| s + e) / n).sqrt();

            let cv = rational::Ratio::from_float(std / mean).unwrap();

            let mut p = pmf.entry(cv).or_insert(0.0);
            *p = logprobs::add(*p, prob);
        }
    }

    let pmf = pmf.iter().map(|(cv, p)| {
        model::pmf::Entry {
            value: cv.numer().to_f64().unwrap().log2() - cv.denom().to_f64().unwrap().log2(),
            prob: *p
        }
    }).collect_vec();

    model::diffexp::PMF::new(pmf)
}
