use itertools::Itertools;
use ndarray::prelude::*;
use statrs::distribution::Multinomial;

use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codebook, Codeword, FeatureID};


pub struct FeatureModel {
    feature_id: u32,
    event_probs: Vec<f64>,
    neighbors: Vec<FeatureID>,
    min_dist: u8
}


impl FeatureModel {
    pub fn new(feature: &str, codebook: &Codebook, xi: &Xi) -> Self {
        assert!(
            codebook.min_dist == 2 || codebook.min_dist == 4,
            "unsupported hamming distance: {}, supported: 2, 4",
            codebook.min_dist
        );

        let feature = codebook.get_id(feature);
        let feature_record = codebook.get_record(feature);

        let xi_call = xi.prob_error(feature_record, None);

        let prob_call_exact = xi_call[(0, 0)];
        let prob_call_mismatch = if codebook.min_dist == 4 {
            xi_call[(1, 0)].ln_add_exp(xi_call[(0, 1)])
        } else {
            LogProb::ln_zero()
        };

        let neighbors_4 = codebook.neighbors(feature, 4);
        let neighbors_2 = if codebook.min_dist == 2 {
            codebook.neighbors(feature, 2)
        } else {
            vec![]
        };

        let mut neighbors = neighbors_4.clone();
        neighbors.extend(&neighbors_2);

        let xi_miscall_4 = codebook.neighbors(feature, 4).map(|n| {
            xi.prob_error(feature_record, Some(codebook.get_record(n)))
        }).collect_vec();
        let xi_miscall_2 = codebook.neighbors(feature, 2).map(|n| {
            xi.prob_error(feature_record, Some(codebook.get_record(n)))
        }).collect_vec();

        let prob_miscall_exact = xi_miscall_4.map(|xi| {
            xi[(2, 2)]
        }).chain(xi_miscall_2.map(|xi| xi[(1, 1)])).collect_vec();

        let prob_miscall_mismatch = if codebook.min_dist == 4 {
            xi_miscall_4.map(|xi| {
                LogProb::ln_sum_exp(&[
                    xi[(2, 1)], xi[(1, 2)], xi[(2, 3)], xi[(3, 2)]
                ])
            });
        } else {
            vec![];
        };

        let prob_total = LogProb::ln_sum_exp(&[
            prob_call_exact,
            prob_call_mismatch,
            LogProb::ln_sum_exp(&prob_miscall_exact),
            LogProb::ln_sum_exp(&prob_miscall_mismatch)
        ]);

        let mut event_probs = vec![
            *Prob::from(prob_call_exact),
            *Prob::from(prob_call_mismatch),
            *Prob::from(prob_total.ln_one_minus_exp()),
        ];
        event_probs.extend(prob_miscall_exact.iter().map(|p| *Prob::from(p)));
        event_probs.extend(prob_miscall_mismatch.iter().map(|p| *Prob::from(p)));

        FeatureModel {
            feature_id: feature,
            event_probs: event_probs,
            neighbors: neighbors,
            min_dist: codebook.min_dist
        }
    }

    pub fn likelihood(
        &self,
        x: u32,
        count: u32,
        count_exact: u32,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> LogProb {
        assert!(count >= count_exact, "bug: count < count_exact");

        let count_mismatch = count - count_exact;
        let count_exact_miscalled = miscalls_exact.slice(s![self.feature_id, ..]).scalar_sum();
        let count_mismatch_miscalled = miscalls_mismatch.slice(
            s![self.feature_id, ..]
        ).scalar_sum();

        if count_exact < count_exact_miscalled {
            return LogProb::ln_zero();
        }
        if count_mismatch < count_mismatch_miscalled {
            return LogProb::ln_zero();
        }
        if (count_exact - count_exact_miscalled + count_mismatch - count_mismatch_miscalled) > x {
            return LogProb::ln_zero();
        }

        let mut event_counts = vec![
            // exact calls
            (count_exact - count_exact_miscalled) as u64,
            // mismatch calls
            (count_mismatch - count_mismatch_miscalled) as u64,
            x.saturating_sub(count_exact + count_mismatch),
        ];
        event_counts.extend(self.neighbors.map(|n| miscalls_exact[(self.feature_id, n)]));
        if self.min_dist == 4 {
            event_counts.extend(self.neighbors.map(|n| miscalls_mismatch[(self.feature_id, n)]));
        }
        if event_counts.iter().sum() != x {
            return LogProb::ln_zero();
        }

        let multinomial = Multinomial::new(self.event_probs, x);

        LogProb::from(Prob(multinomial.pmf(&event_counts)))
    }
}


type Miscalls = Array2<u32>;


/// Basic model for the probability of making i 1-0 and j 0-1 errors.
struct Xi {
    p0: Vec<LogProb>,
    p1: Vec<LogProb>,
    max_err: u8
}


impl Xi {
    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `p0` - error probabilities for 0-1 errors
    /// * `p1` - error probabilities for 1-0 errors
    /// * `max_err` - maximum number of errors considered
    fn new(p0: &[Prob], p1: &[Prob], max_err: u8) -> Self {
        let tolog = |&p| LogProb::from(p);
        Xi {
            p0: p0.iter().map(&tolog).collect_vec(),
            p1: p1.iter().map(&tolog).collect_vec(),
            max_err: max_err
        }
    }

    fn prob_error(&self &codebook::Record, target: Option<&codebook::Record>) -> Array2<LogProb> {
        let mask = if let Some(target) = target {
            Some(source.diff(target))
        } else {
            None
        };

        self.calc(source.codeword(), mask.as_ref())
    }

    fn calc(&self, codeword: &Codeword, mask: Option<&BitVec>) -> Array2<LogProb> {
        let n = self.max_err as usize;
        let mut curr = Array::from_elem((n + 1, n + 1), LogProb::ln_zero());
        let mut prev = Array::from_elem((n + 1, n + 1), LogProb::ln_zero());
        prev[(1, 1)] = LogProb::ln_one();

        for k in 0..codeword.len() {
            for i in 1..curr.shape()[0] {
                for j in 1..curr.shape()[1] {

                    let (i_err, j_err, p_err) = if codeword.get(k).unwrap() {
                        (i - 1, j, self.p1[k])
                    } else {
                        (i, j - 1, self.p0[k])
                    };

                    let mut p = p_err.ln_one_minus_exp() + prev[(i, j)];

                    if mask.map_or(true, |mask| mask.get(k).unwrap()) {
                        p = p.ln_add_exp(p_err + prev[(i_err, j_err)]);
                    }

                    curr[(i, j)] = p;
                }
            }
            mem::swap(&mut prev, &mut curr);
        }

        prev.slice(s![1.., 1..]).to_owned()
    }
}
