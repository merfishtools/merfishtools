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
    for ((a_mean, a_prob), (b_mean, b_prob)) in a.iter().cartesian_product(b.iter()) {
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

    PMF::new(pmf.iter().map(|(fc, prob)| (*fc.numer() as f64 / *fc.denom() as f64, *prob)).collect_vec())
}


pub fn differential_expression_pep(pmf: &PMF, min_fc: LogFC) -> LogProb {
    let probs = pmf.iter().filter(|&&(fc, _)| fc >= min_fc).map(|&(_, prob)| prob).collect_vec();
    logprobs::ln_1m_exp(logprobs::log_prob_sum(&probs))
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
