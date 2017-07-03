use itertools::Itertools;
use ordered_float::NotNaN;

use bio::stats::LogProb;
use bio::stats::probs;

use model;


pub type CDF = probs::cdf::CDF<NotNaN<f64>>;


/// Calculate CDF of expression.
///
/// Returns a tuple of CDF and the naive estimate.
pub fn cdf(feature: &str, count: u32, count_exact: u32, model: &Box<model::readout::Model>, window_width: u32) -> (CDF, u32) {
    let readout_model = model::Readout::new(feature, model, window_width);
    let (xmin, xmax) = readout_model.window(count, count_exact);
    let likelihoods = (xmin..xmax + 1).map(|x| {
        readout_model.likelihood(x, count, count_exact)
    }).collect_vec();
    // calculate (marginal / flat_prior)
    let marginal = LogProb::ln_sum_exp(&likelihoods);

    (CDF::from_pmf(
        likelihoods.iter().enumerate().filter_map(|(x, lh)| {
            let prob = lh - marginal;
            if prob >= model::MIN_PROB {
                Some(probs::cdf::Entry {
                    value: NotNaN::new((xmin + x as u32) as f64).unwrap(),
                    prob: prob
                })
            }
            else {
                None
            }
        }).collect_vec()
    ), readout_model.naive_estimate(count))
}



#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use nalgebra::ApproxEq;

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

    fn setup_mhd2() -> Box<model::readout::Model> {
        model::readout::new_model(
            &[Prob(0.04); 14],
            &[Prob(0.1); 14],
            io::codebook::Codebook::from_file("test/codebook/simulated-MHD2.txt").unwrap()
        )
    }

    #[test]
    fn test_cdf() {
        let readout = setup();
        let (cdf, _) = cdf(GENE, 25, 10, &readout, 100);

        let total = cdf.total_prob();
        println!("{:?}", cdf);
        println!("{}", *total);
        assert!(total.approx_eq(&-0.0000011368907495423741));
    }

    #[test]
    fn test_cdf2() {
        let readout = setup();
        let (cdf, _) = cdf(GENE, 176, 25, &readout, 100);

        let total = cdf.total_prob();
        println!("{:?}", cdf);
        println!("{}", *total);
        assert!(total.approx_eq(&-0.0000035876739698048574));
    }

    #[test]
    fn test_cdf3() {
        let calls = 40;
        let miscalls = 8;
        let count = calls + miscalls;

        let readout = setup_mhd2();
        let (cdf, _) = cdf("COL7A1", count, count, &readout, 100);

        println!("{} {:?}", cdf.map().unwrap(), cdf);
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
