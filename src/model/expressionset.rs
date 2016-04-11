use model;

pub type MeanExpression = f64;
pub type CDF = model::dist::CDF<MeanExpression>;


pub fn cdf(expression_pmfs: &[model::expression::CDF], pseudocounts: f64) -> CDF {
    let test = model::meanvar::cdf(expression_pmfs, |mean, var| (mean, var));
    println!("DEBUG {}", test.total_prob());
    model::meanvar::cdf(expression_pmfs, |mean, _| mean + pseudocounts)
}


#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use itertools::Itertools;
    use nalgebra::ApproxEq;
    use bio::stats::logprobs;
    use bio::stats::logprobs::Prob;

    use super::*;
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
    fn test_cdf() {
        let readout = setup();
        let cdfs = [
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        println!("{:?}", cdfs[0]);
        let cdf = cdf(&cdfs, 0.0);
        println!("{:?}", cdf);

        let total = cdf.total_prob();

        assert_relative_eq!(total, 0.0, epsilon = 0.0002);
    }
}
