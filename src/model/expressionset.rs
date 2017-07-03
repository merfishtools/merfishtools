use ordered_float::NotNaN;

use bio::stats::probs;

use model;

pub type MeanExpression = NotNaN<f64>;
pub type CDF = probs::cdf::CDF<MeanExpression>;


pub fn cdf(expression_cdfs: &[model::expression::CDF], pseudocounts: f64) -> CDF {
    assert!(pseudocounts != 0.0, "pseudocounts may not be zero");
    let pseudocounts = NotNaN::new(pseudocounts).unwrap();
    let cdf = model::meanvar::cdf(expression_cdfs, |mean, _| mean + pseudocounts);
    cdf.reduce().sample(1000)
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use bio::stats::{Prob, LogProb};

    use super::*;
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
    fn test_cdf() {
        let readout = setup();
        let cdfs = [
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            model::expression::cdf(GENE, 5, 5, &readout, 100).0
        ];
        println!("{:?}", cdfs[0]);
        let cdf = cdf(&cdfs, 0.0);
        println!("{:?}", cdf);

        let total = cdf.total_prob();

        assert_relative_eq!(*total, *LogProb::ln_one(), epsilon = 0.0002);
    }
}
