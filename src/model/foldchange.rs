use num::rational;
use itertools::Itertools;

use model;


pub type LogFC = f64;
pub type FC = rational::Ratio<u64>;

pub type CDF = model::dist::CDF<LogFC>;


/// PMF of log2 fold change of a vs b. Specifically, we calculate log2((mean(a) + c) / (mean(b) + c))
pub fn cdf(a: &model::expressionset::CDF, b: &model::expressionset::CDF) -> CDF {
    let mut pmf = Vec::new();
    let a_pmf = a.iter_pmf().collect_vec();
    let b_pmf = b.iter_pmf().collect_vec();

    for (a_mean, a_prob) in a_pmf {
        for &(b_mean, b_prob) in b_pmf.iter() {
            let log2fc = (a_mean).log2() - (b_mean).log2();
            pmf.push((log2fc, a_prob + b_prob));
        }
    }
    model::dist::CDF::from_pmf(pmf).reduce().sample(100)
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]
    use super::*;

    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("test/codebook/140genesData.1.txt", 4).unwrap().codebook())
    }

    #[test]
    fn test_foldchange_pmf() {
        let readout = setup();
        let cdfs1 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        let cdfs2 = [
            model::expression::cdf(GENE, 50, 50, &readout, 100),
            model::expression::cdf(GENE, 50, 50, &readout, 100),
            model::expression::cdf(GENE, 50, 50, &readout, 100),
            model::expression::cdf(GENE, 50, 50, &readout, 100)
        ];
        let cdf1 = model::expressionset::cdf(&cdfs1, 0.0);
        let cdf2 = model::expressionset::cdf(&cdfs2, 0.0);

        let cdf = cdf(&cdf2, &cdf1);

        let total = cdf.total_prob();
        let fc = 2.0f64.powf(cdf.expected_value());

        println!("ev={}", fc);
        assert_relative_eq!(total, 0.0, epsilon = 0.0002);
        assert!(fc >= 9.0);
    }
}
