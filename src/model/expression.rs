use itertools::Itertools;

use bio::stats::logprobs;

use model;


pub type PMF = model::pmf::PMF<f32>;


pub fn pmf(count: u32, count_exact: u32, readout_model: &model::Readout) -> PMF {
    let (xmin, xmax) = readout_model.window(count);
    let likelihoods = (xmin..xmax + 1).map(|x| {
        readout_model.likelihood(x, count, count_exact)
    }).collect_vec();
    // calculate (marginal / flat_prior)
    let marginal = logprobs::log_prob_sum(&likelihoods);

    // TODO trim
    PMF::new(
        likelihoods.iter().enumerate().map(|(x, lh)| {
            model::pmf::Entry{
                value: (xmin + x as u32) as f32,
                prob: lh - marginal
            }
        }).filter(|e| e.prob >= model::MIN_PROB).collect_vec()
    )
}



#[cfg(test)]
mod tests {
    #![allow(non_upper_case_globals)]

    use bio::stats::logprobs::Prob;

    use super::*;
    use model;


    const N: u8 = 16;
    const m: u8 = 4;
    const p0: Prob = 0.04;
    const p1: Prob = 0.1;
    const dropout_rate: Prob = 0.15;

    fn setup() -> model::Readout {
        model::Readout::new(N, m, p0, p1, dropout_rate)
    }

    #[test]
    fn test_pmf() {
        let readout = setup();
        let pmf = pmf(0, 0, &readout);

        println!("{:?}", pmf);
        assert_eq!(pmf.map(), 0.0);
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
