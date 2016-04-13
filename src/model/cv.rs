use model;


pub type CV = f64;
pub type CDF = model::dist::CDF<CV>;


pub fn cdf(cdfs: &[model::expressionset::CDF]) -> CDF {
    model::meanvar::cdf(cdfs, |mean, var| var.sqrt() / mean).reduce().sample(100)
}


#[cfg(test)]
mod tests {
    use super::*;

    use itertools::Itertools;
    use bio::stats::logprobs;

    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("evaluation/codebook/140genesData.1.txt", 4).unwrap().codebook())
    }

    #[test]
    fn test_cdf() {
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
        println!("{} {}", cdfs1[0].expected_value(), cdfs2[0].expected_value());

        let cdf = cdf(&[cdf1, cdf2]);

        let total = cdf.total_prob();
        println!("{:?}", cdf);

        assert!(total <= 0.0);
        assert_relative_eq!(total, 0.0, epsilon = 0.0002);
        assert_relative_eq!(cdf.expected_value(), 1.14, epsilon = 0.02);
    }
}
