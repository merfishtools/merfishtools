// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use ordered_float::NotNan;

use bio::stats::probs;
use bio::stats::LogProb;

use crate::model;

pub type DiffexpMeasure = NotNan<f64>;
pub type CDF = probs::cdf::CDF<DiffexpMeasure>;

/// An estimate of differential expression.
#[derive(Debug)]
pub struct Estimate {
    pub differential_expression_pep: LogProb,
    pub differential_expression_bf: f64,
    pub map: DiffexpMeasure,
    pub credible_interval: (DiffexpMeasure, DiffexpMeasure),
}

/// Posterior error probability for differential expression.
pub fn pep(cdf: &CDF, max_null_value: DiffexpMeasure) -> LogProb {
    cdf.get(&max_null_value).unwrap()
}

/// Bayes factor for differential expression.
pub fn bayes_factor(cdf: &CDF, max_null_value: DiffexpMeasure) -> model::bayes::BayesFactor {
    let m0 = cdf.get(&max_null_value).unwrap();
    let m1 = m0.ln_one_minus_exp();

    (m1 - m0).exp()
}

pub fn estimate(cdf: &CDF, max_fc: DiffexpMeasure) -> Estimate {
    let ci = cdf.credible_interval(0.95).expect("bug: empty CDF");
    Estimate {
        differential_expression_pep: pep(cdf, max_fc),
        differential_expression_bf: bayes_factor(cdf, max_fc),
        map: *cdf.map().expect("bug: empty CDF"),
        credible_interval: (*ci.start, *ci.end),
    }
}
