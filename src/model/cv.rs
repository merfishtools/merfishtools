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

    for groups in combinations {
        let prob = groups.iter().map(|group| group.prob).fold(0.0, |s, p| s + p);
        if prob >= model::MIN_PROB {
            // mean
            let mean = {
                let summands = groups.iter().map(|group| {
                    group.prob +
                    (*group.value.numer() as f64 / *group.value.denom() as f64).ln()
                }).collect_vec();
                logprobs::sum(&summands).exp()
            };
            // standard deviation
            let std = {
                let summands = groups.iter().map(|group| {
                    group.prob +
                    (*group.value.numer() as f64 / *group.value.denom() as f64 - mean).ln().powi(2)
                }).collect_vec();
                logprobs::sum(&summands).exp().sqrt()
            };
            let cv = rational::Ratio::from_float(std / mean).unwrap();

            let mut p = pmf.entry(cv).or_insert(0.0);
            *p = logprobs::add(*p, prob);

            //pmf.push(model::pmf::Entry { value: cv, prob: prob });
        }
    }

    let pmf = pmf.iter().map(|(cv, p)| {
        model::pmf::Entry { value: cv.numer().to_f64().unwrap() / cv.denom().to_f64().unwrap(), prob: *p }
    }).collect_vec();

    model::diffexp::PMF::new(pmf)
}
