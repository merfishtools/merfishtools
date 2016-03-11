use std::collections;

use num::rational;
use itertools::Itertools;
use bio::stats::logprobs;
use bio::stats::logprobs::LogProb;

use model;


pub type LogFC = f64;
pub type FC = rational::Ratio<u64>;
pub type PMF = model::pmf::PMF<LogFC>;


pub fn pmf(a: &model::expressionset::PMF, b: &model::expressionset::PMF) -> PMF {
    let mut pmf = collections::HashMap::new();
    for (a, b) in a.iter().cartesian_product(b.iter()) {
        // add pseudocount
        let fc = (b.value + rational::Ratio::from_integer(1)) / (a.value + rational::Ratio::from_integer(1));
        let posterior_prob = b.prob + a.prob;

        if pmf.contains_key(&fc) {
            let p = pmf.get_mut(&fc).unwrap();
            *p = logprobs::add(*p, posterior_prob);
        }
        else {
            pmf.insert(fc, posterior_prob);
        }
    }

    let mut pmf = pmf.iter().map(|(fc, prob)| {
        model::pmf::Entry { value: (*fc.numer() as f64 + 1.0).log2() - (*fc.denom() as f64 + 1.0).log2(), prob: *prob }
    }).collect_vec();
    pmf.sort_by(|a, b| a.value.partial_cmp(&b.value).unwrap());

    PMF::new(pmf)
}


impl PMF {
    pub fn differential_expression_pep(&self, min_fc: LogFC) -> LogProb {
        let probs = self.iter().filter(|e| e.value.abs() > min_fc).map(|e| e.prob).collect_vec();
        logprobs::ln_1m_exp(logprobs::sum(&probs))
    }
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use itertools::Itertools;
    use nalgebra::ApproxEq;
    use bio::stats::logprobs;
    use bio::stats::logprobs::Prob;

    use model;
    use io;

    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;
    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("evaluation/codebook/140genesData.1.txt", 4).unwrap().codebook())
    }

    #[test]
    fn test_pmf() {
        let readout = setup();
        let pmfs1 = [
            model::expression::pmf(GENE, 5, 5, &readout),
            model::expression::pmf(GENE, 5, 5, &readout),
            model::expression::pmf(GENE, 5, 5, &readout),
            model::expression::pmf(GENE, 5, 5, &readout)
        ];
        let pmfs2 = [
            model::expression::pmf(GENE, 50, 50, &readout),
            model::expression::pmf(GENE, 50, 50, &readout),
            model::expression::pmf(GENE, 50, 50, &readout),
            model::expression::pmf(GENE, 50, 50, &readout)
        ];
        let pmf1 = model::expressionset::pmf(&pmfs1);
        let pmf2 = model::expressionset::pmf(&pmfs2);

        let pmf = pmf(&pmf1, &pmf2);

        let total = logprobs::sum(&pmf.iter().map(|fc| fc.prob).collect_vec());
        let fc = 2.0f64.powf(pmf.expected_value());

        println!("{:?}", total);
        println!("ev={}", fc);
        assert!(total.approx_eq(&-0.000014479671117229032));
    }
}
