// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use model;


pub type CV = f64;
pub type CDF = model::dist::CDF<CV>;


/// Calculate CDF for differential expression.
pub fn cdf(cdfs: &[model::expressionset::CDF]) -> CDF {
    model::meanvar::cdf(cdfs, |mean, var| var.sqrt() / mean).reduce().sample(100)
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("test/codebook/140genesData.1.txt", 4).unwrap().codebook())
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

    #[test]
    fn test_zero_cv() {
        let readout = setup();
        let cdfs1 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        let cdfs2 = [
            model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100),
            //model::expression::cdf(GENE, 5, 5, &readout, 100)
        ];
        let cdf1 = model::expressionset::cdf(&cdfs1, 0.0);
        let cdf2 = model::expressionset::cdf(&cdfs2, 0.0);
        println!("{} {}", cdfs1[0].expected_value(), cdfs2[0].expected_value());
        println!("{:?}", cdf1.iter_pmf().collect_vec());

        let cdf = cdf(&[cdf1, cdf2]);

        assert_relative_eq!(*cdf.map(), 0.0);
    }

    #[test]
    fn test_ahnak2() {
        // let mut cdfs1 = vec![
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 2, 2, &readout, 100),
        //     model::expression::cdf(gene, 3, 3, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 2, 2, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 6, 6, &readout, 100),
        //     model::expression::cdf(gene, 2, 2, &readout, 100),
        //     model::expression::cdf(gene, 4, 4, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        // ];
        // cdfs1.
        // let mut cdfs2 = vec![
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 1, 1, &readout, 100),
        //     model::expression::cdf(gene, 2, 2, &readout, 100)
        // ];
        let cdfs1 = io::cdf::expression::Reader::from_file("test/expression_cdf/1.txt").unwrap().cdfs();
        let cdfs2 = io::cdf::expression::Reader::from_file("test/expression_cdf/2.txt").unwrap().cdfs();
        let cdfs3 = io::cdf::expression::Reader::from_file("test/expression_cdf/3.txt").unwrap().cdfs();

        let cdfs1 = cdfs1.get("AHNAK2").unwrap();
        let cdfs2 = cdfs2.get("AHNAK2").unwrap();
        let cdfs3 = cdfs3.get("AHNAK2").unwrap();

        let cdf1 = model::expressionset::cdf(&cdfs1, 1.0);
        let cdf2 = model::expressionset::cdf(&cdfs2, 1.0);
        let cdf3 = model::expressionset::cdf(&cdfs3, 1.0);
        println!("{} {} {}", cdf1.map(), cdf2.map(), cdf3.map());
        println!("{} {} {}", cdf1.expected_value(), cdf2.expected_value(), cdf3.expected_value());

        let cdf = cdf(&[cdf1, cdf2]);
        println!("{}", cdf.map());
        println!("{}", cdf.expected_value());

        println!("{:?}", cdf.estimate(0.5));

        assert!(false);
    }
}
