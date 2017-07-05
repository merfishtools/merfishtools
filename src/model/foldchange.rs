use ordered_float::NotNaN;
use num::rational;
use itertools::Itertools;

use bio::stats::probs;

use model;


pub type LogFC = NotNaN<f64>;
pub type FC = rational::Ratio<u64>;

pub type CDF = probs::cdf::CDF<LogFC>;


/// PMF of log2 fold change of a vs b. Specifically, we calculate log2((mean(a) + c) / (mean(b) + c))
pub fn cdf(a: &model::expressionset::CDF, b: &model::expressionset::CDF) -> CDF {
    let mut pmf = Vec::new();
    let a_pmf = a.iter_pmf().collect_vec();
    let b_pmf = b.iter_pmf().collect_vec();

    for a in a_pmf {
        let a_mean = *a.value;
        for b in b_pmf.iter() {
            let b_mean = *b.value;
            // the PMFs should not contain NaNs.
            let log2fc = NotNaN::new((a_mean).log2() - (b_mean).log2()).unwrap();
            pmf.push(probs::cdf::Entry { value: log2fc, prob: a.prob + b.prob });
        }
    }
    CDF::from_pmf(pmf).reduce().sample(100)
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use bio::stats::{Prob, LogProb};

    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(
            &[Prob(0.04); 16],
            &[Prob(0.1); 16],
            io::codebook::Codebook::from_file("test/codebook/140genesData.1.txt").unwrap()
        )
    }

    #[test]
    fn test_foldchange_pmf() {
        let readout = setup();
        let cdfs1 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0
        ];
        let cdfs2 = [
            model::expression::cdf(GENE, 50, 50, &readout, 100).0,
            model::expression::cdf(GENE, 50, 50, &readout, 100).0,
            model::expression::cdf(GENE, 50, 50, &readout, 100).0,
            model::expression::cdf(GENE, 50, 50, &readout, 100).0
        ];
        let cdf1 = model::expressionset::cdf(&cdfs1, 0.000000001);
        let cdf2 = model::expressionset::cdf(&cdfs2, 0.000000001);

        let cdf = cdf(&cdf2, &cdf1);

        let total = cdf.total_prob();
        let fc = 2.0f64.powf(**cdf.map().unwrap());

        println!("map={}", fc);
        assert_relative_eq!(*total, *LogProb::ln_one(), epsilon = 0.002);
        assert!(fc >= 9.0);
    }
}
