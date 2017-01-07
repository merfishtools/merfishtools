use itertools::Itertools;

use bio::stats::logprobs;

use model;


pub type CDF = model::dist::CDF<f64>;


pub fn cdf(feature: &str, count: u32, count_exact: u32, model: &Box<model::readout::Model>, window_width: u32) -> CDF {
    let readout_model = model::Readout::new(feature, model, window_width);
    let (xmin, xmax) = readout_model.window(count, count_exact);
    let likelihoods = (xmin..xmax + 1).map(|x| {
        readout_model.likelihood(x, count, count_exact)
    }).collect_vec();
    // calculate (marginal / flat_prior)
    let marginal = logprobs::sum(&likelihoods);

    model::dist::CDF::from_pmf(
        likelihoods.iter().enumerate().filter_map(|(x, lh)| {
            let prob = lh - marginal;
            if prob >= model::MIN_PROB {
                Some(((xmin + x as u32) as f64, prob))
            }
            else {
                None
            }
        }).collect_vec()
    )
}



#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use nalgebra::ApproxEq;

    use super::*;
    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("test/codebook/140genesData.1.txt", 4).unwrap().codebook())
    }

    fn setup_mhd2() -> Box<model::readout::Model> {
        model::readout::new_model(14, 4, 0.04, 0.1, 2, io::codebook::Reader::from_file("test/codebook/simulated-MHD2.txt", 2).unwrap().codebook())
    }

    #[test]
    fn test_cdf() {
        let readout = setup();
        let cdf = cdf(GENE, 25, 10, &readout, 100);

        let total = cdf.total_prob();
        println!("{:?}", cdf);
        println!("{}", total);
        assert!(total.approx_eq(&-0.0000011368907495423741));
    }

    #[test]
    fn test_cdf2() {
        let readout = setup();
        let cdf = cdf(GENE, 176, 25, &readout, 100);

        let total = cdf.total_prob();
        println!("{:?}", cdf);
        println!("{}", total);
        assert!(total.approx_eq(&-0.0000035876739698048574));
    }

    #[test]
    fn test_cdf3() {
        let calls = 40;
        let miscalls = 8;
        let count = calls + miscalls;

        let readout = setup_mhd2();
        let cdf = cdf("COL7A1", count, count, &readout, 100);

        println!("{:?} {} {:?}", cdf.expected_value(), cdf.map(), cdf);
    }
}
/*
    #[test]
    fn test_likelihood() {
        let readout = setup();
        assert!(likelihood(0, 0, 0, &readout).exp().approx_eq(&1.0));
        assert!(likelihood(1, 1, 1, &readout).exp().approx_eq(&0.4019988717840602));
    }

    #[test]
    fn test_posterior_prob() {
        let readout = setup();
        let expression = Expression::new(100, 99, &readout);
        println!("{:?}", (90..110).map(|x| expression.posterior_prob(x).exp()).collect_vec());
        assert!(false);
        let expression = Expression::new(5, 5, &readout);
        // check if x=5 yields highest probability
        assert_eq!((0..21).sorted_by(|&x, &y| {
            expression.posterior_prob(x).partial_cmp(&expression.posterior_prob(y)).unwrap()
        })[20], 5);
        // check that expression beyond window yields zero probability
        assert!(expression.posterior_prob(100).exp().approx_eq(&0.0));
    }

    // #[test]
    // fn test_set_expected_value() {
    //     let readout = setup();
    //     let expressions = [
    //         Expression::new(5, 1, &readout),
    //         Expression::new(5, 1, &readout),
    //         Expression::new(5, 1, &readout),
    //         Expression::new(1, 0, &readout)
    //     ];
    //     let expression_set = ExpressionSet::new(&expressions);
    //     // calculate expected value directly
    //     let expected_value = expressions.iter().map(|e| e.expected_value()).fold(0.0, |s, e| s + e) / expressions.len() as f64;
    //     assert!(expression_set.expected_value().approx_eq(&expected_value));
    // }
}
*/
