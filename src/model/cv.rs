// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use ordered_float::NotNaN;

use bio::stats::probs;

use model;


pub type CV = NotNaN<f64>;
pub type CDF = probs::cdf::CDF<CV>;


/// Calculate CDF for differential expression.
pub fn cdf(cdfs: &[model::expressionset::CDF]) -> CDF {
        model::meanvar::cdf(cdfs, |mean, var| NotNaN::new(var.sqrt()).unwrap() / mean).reduce().sample(100)
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

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
    fn test_cdf() {
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
        println!("{} {}", cdfs1[0].map().unwrap(), cdfs2[0].map().unwrap());

        let cdf = cdf(&[cdf1, cdf2]);

        let total = cdf.total_prob();
        println!("{:?}", cdf);

        assert!(*total <= 0.0);
        assert_relative_eq!(*total, 0.0, epsilon = 0.005);
        assert_relative_eq!(**cdf.map().unwrap(), 1.14, epsilon = 0.02);
    }

    #[test]
    fn test_zero_cv() {
        let readout = setup();
        let cdfs1 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        let cdfs2 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100).0,
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        let cdf1 = model::expressionset::cdf(&cdfs1, 0.0);
        let cdf2 = model::expressionset::cdf(&cdfs2, 0.0);
        println!("{:?}", cdf1.iter_pmf().collect_vec());

        let cdf = cdf(&[cdf1, cdf2]);

        assert_relative_eq!(**cdf.map().unwrap(), 0.0);
    }

    #[test]
    fn test_ahnak2() {
        let cdfs1 = io::cdf::expression::Reader::from_file("test/expression_cdf/1.txt").unwrap().cdfs();
        let cdfs2 = io::cdf::expression::Reader::from_file("test/expression_cdf/2.txt").unwrap().cdfs();
        let cdfs3 = io::cdf::expression::Reader::from_file("test/expression_cdf/3.txt").unwrap().cdfs();

        let cdfs1 = cdfs1.get("AHNAK2").unwrap();
        let cdfs2 = cdfs2.get("AHNAK2").unwrap();
        let cdfs3 = cdfs3.get("AHNAK2").unwrap();

        let cdf1 = model::expressionset::cdf(&cdfs1, 1.0);
        let cdf2 = model::expressionset::cdf(&cdfs2, 1.0);
        let cdf3 = model::expressionset::cdf(&cdfs3, 1.0);
        println!("{} {} {}", cdf1.map().unwrap(), cdf2.map().unwrap(), cdf3.map().unwrap());

        let cdf = cdf(&[cdf1, cdf2, cdf3]);
        println!("{}", cdf.map().unwrap());

        println!("{:?}", model::diffexp::estimate(&cdf, NotNaN::new(0.5).unwrap()));

        // TODO add proper test case
    }
}
