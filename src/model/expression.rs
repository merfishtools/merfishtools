use itertools::Itertools;
use ordered_float::NotNaN;

use bio::stats::LogProb;
use bio::stats::probs;

use model;
use io::codebook::FeatureID;


pub type CDF = probs::cdf::CDF<u32>;
pub type NormalizedCDF = probs::cdf::CDF<NotNaN<f64>>;


/// Calculate CDF of expression.
///
/// Returns a tuple of CDF and the MAP estimate.
pub fn cdf(feature: FeatureID, model: &mut model::readout::JointModel) -> (CDF, u32) {
    let (xmin, xmax) = model.window(feature);
    let likelihoods = (xmin..xmax + 1).map(|x| {
        model.likelihood(feature, x)
    }).collect_vec();
    // calculate (marginal / flat_prior)
    let marginal = LogProb::ln_sum_exp(&likelihoods);

    let cdf = CDF::from_pmf(
        likelihoods.iter().enumerate().filter_map(|(x, lh)| {
            let prob = lh - marginal;
            if prob >= model::MIN_PROB {
                Some(probs::cdf::Entry {
                    value: xmin + x as u32,
                    prob
                })
            }
            else {
                None
            }
        }).collect_vec()
    );

    let map = model.map_estimate(feature);
    let cdf_map = *cdf.map().expect(&format!("bug: empty CDF, xmin={}, xmax={}", xmin, xmax));

    if cdf_map != map {
        debug!("EM-MAP and CDF MAP do not agree: {}!={}", map, cdf_map);
    }

    // assert_eq!(
    //     cdf_map, map,
    //     "bug: CDF-derived MAP and EM MAP are not the same: {} != {}", cdf_map, map
    // );

    (cdf, map)
}



#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

//    use nalgebra::ApproxEq;
//
//    use bio::stats::{Prob, LogProb};
//
//    use super::*;
//    use model;
//    use io;

/*

    const GENE: &'static str = "COL5A1";

    fn setup() -> model::readout::JointModel {
        model::readout::new_model(
            &[Prob(0.04); 16],
            &[Prob(0.1); 16],
            io::codebook::Codebook::from_file("tests/codebook/140genesData.1.txt").unwrap()
        )
    }

    fn setup_mhd2() -> Box<model::readout::Model> {
        model::readout::new_model(
            &[Prob(0.04); 14],
            &[Prob(0.1); 14],
            io::codebook::Codebook::from_file("tests/codebook/simulated-MHD2.txt").unwrap()
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
    }*/
}
