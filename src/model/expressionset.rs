use model;

pub type MeanExpression = f64;
pub type CDF = model::dist::CDF<MeanExpression>;


pub fn cdf(expression_cdfs: &[model::expression::CDF], pseudocounts: f64) -> CDF {
    model::meanvar::cdf(expression_cdfs, |mean, _| mean + pseudocounts).reduce().sample(1000)
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
