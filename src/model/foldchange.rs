use std::collections;

use num::rational;
use itertools::Itertools;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


pub type LogFC = f64;
pub type FC = rational::Ratio<u32>;
pub type PMF = model::pmf::PMF<LogFC>;


pub fn pmf(a: &model::expressionset::PMF, b: &model::expressionset::PMF) -> PMF {
    let mut pmf = collections::HashMap::new();
    for (&(a_mean, a_prob), &(b_mean, b_prob)) in a.iter().cartesian_product(b.iter()) {
        // add pseudocount
        let fc = (b_mean + rational::Ratio::from_integer(1)) / (a_mean + rational::Ratio::from_integer(1));
        let posterior_prob = b_prob + a_prob;

        if pmf.contains_key(&fc) {
            let p = pmf.get_mut(&fc).unwrap();
            *p = logprobs::log_prob_add(*p, posterior_prob);
        }
        else {
            pmf.insert(fc, posterior_prob);
        }
    }

    PMF::new(pmf.iter().map(|(fc, prob)| (*fc.numer() as f64 / *fc.denom() as f64, *prob)).collect_vec())
}


impl PMF {
    pub fn differential_expression_pep(&self, min_fc: LogFC) -> LogProb {
        let probs = self.iter().filter(|&&(fc, _)| fc >= min_fc).map(|&(_, prob)| prob).collect_vec();
        logprobs::ln_1m_exp(logprobs::log_prob_sum(&probs))
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use itertools::Itertools;
    use nalgebra::ApproxEq;
    use bio::stats::logprobs::{Prob, log_prob_sum};

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
        let pmfs1 = [
            model::expression::pmf(5, 1, &readout),
            model::expression::pmf(10, 1, &readout),
            model::expression::pmf(3, 1, &readout),
            model::expression::pmf(24, 1, &readout)
        ];
        let pmfs2 = [
            model::expression::pmf(50, 1, &readout),
            model::expression::pmf(100, 1, &readout),
            model::expression::pmf(30, 1, &readout),
            model::expression::pmf(240, 1, &readout)
        ];
        let scales = [1.0, 1.0, 1.0, 1.0];
        let pmf1 = model::expressionset::pmf(&pmfs1, &scales);
        let pmf2 = model::expressionset::pmf(&pmfs2, &scales);

        let pmf = pmf(&pmf1, &pmf2);


        let total = log_prob_sum(&pmf.iter().map(|&(_, prob)| prob).collect_vec());

        println!("{:?}", total);
        println!("ev={}", pmf.expected_value());
        assert!(total.approx_eq(&-0.00000948431106451153));
    }
}
