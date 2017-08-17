use std::mem;
use std::cell::RefCell;
use std::cmp;

use rand;
use rand::Rng;
use itertools::Itertools;
use ndarray::prelude::*;
use ndarray;
use bit_vec::BitVec;
use rgsl::randist::multinomial::multinomial_pdf;

use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codebook, Codeword, FeatureID};

use model::readout::{Expressions, Counts, Miscalls, Xi};


#[derive(Debug)]
pub enum FeatureType {
    Expressed,
    NotExpressed,
    Noise
}


impl FeatureType {
    pub fn is_expressed(&self) -> bool {
        match self {
            &FeatureType::NotExpressed => false,
            _ => true
        }
    }
}


pub struct FeatureModel {
    pub(crate) feature_id: FeatureID,
    event_probs: Vec<f64>,
    neighbors: Vec<FeatureID>,
    min_dist: u8,
    pub(crate) counts: Counts,
    event_counts: RefCell<Vec<u32>>,
    pub(crate) feature_type: FeatureType
}


impl FeatureModel {
    pub fn noise(feature_id: FeatureID, codebook: &Codebook, xi: &Xi) -> Self {
        assert!(
            codebook.min_dist == 2 || codebook.min_dist == 4,
            "unsupported hamming distance: {}, supported: 2, 4",
            codebook.min_dist
        );




        // miscalls generated from noise
        let xi_miscall = codebook.records().iter().map(
            |feat_rec| xi.prob_error(codebook.noise(), Some(feat_rec), codebook.m as usize + 1)
        ).collect_vec();
        let neighbors = codebook.records().iter().map(
            |rec| codebook.get_id(rec.name())
        ).collect_vec();

        let mut miscall_probs = xi_miscall.iter().map(|xi| {
            match codebook.m {
                4 => xi[(0, 4)],
                8 => xi[(0, 8)],
                6 => xi[(0, 6)],
                _ => panic!("unsupported number of 1-bits")
            }
        }).collect_vec();
        if codebook.min_dist == 4 {
            miscall_probs.extend(xi_miscall.iter().map(|xi| {
                match codebook.m {
                    4 => LogProb::ln_sum_exp(&[xi[(0, 3)], xi[(0, 5)]]),
                    6 => LogProb::ln_sum_exp(&[xi[(0, 5)], xi[(0, 7)]]),
                    8 => LogProb::ln_sum_exp(&[xi[(0, 7)], xi[(0, 9)]]),
                    _ => panic!("unsupported number of 1-bits")
                }
            }));
        }

        let total_miscall_prob = LogProb::ln_sum_exp(&miscall_probs);

        // no call exact or call mismatch events for noise
        let mut event_probs = vec![0.0; 2];
        // dropout event means noise that is not generating miscalls
        event_probs.push(*Prob::from(total_miscall_prob.ln_one_minus_exp()));
        // miscalls
        event_probs.extend(miscall_probs.into_iter().map(|p| *Prob::from(p)));

        FeatureModel {
            feature_id: feature_id,
            event_probs: event_probs,
            neighbors: neighbors,
            min_dist: codebook.min_dist,
            counts: Counts { exact: 0, mismatch: 0 },
            event_counts: RefCell::new(Vec::new()),
            feature_type: FeatureType::Noise
        }
    }

    pub fn new(
        feature: &str,
        counts: Counts,
        codebook: &Codebook,
        xi: &Xi
    ) -> Self {
        assert!(
            codebook.min_dist == 2 || codebook.min_dist == 4,
            "unsupported hamming distance: {}, supported: 2, 4",
            codebook.min_dist
        );

        let feature = codebook.get_id(feature);
        let feature_record = codebook.record(feature);

        let xi_call = xi.prob_error(feature_record, None, 1);

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

        let xi_miscall_4 = codebook.neighbors(feature, 4).iter().map(|&n| {
            xi.prob_error(feature_record, Some(codebook.record(n)), 3)
        }).collect_vec();
        let xi_miscall_2 = if codebook.min_dist == 2 {
            codebook.neighbors(feature, 2).iter().map(|&n| {
                xi.prob_error(feature_record, Some(codebook.record(n)), 1)
            }).collect_vec()
        } else {
            vec![]
        };

        let prob_miscall_exact = xi_miscall_4.iter().map(|xi| {
            xi[(2, 2)]
        }).chain(xi_miscall_2.iter().map(|xi| xi[(1, 1)])).collect_vec();

        let prob_miscall_mismatch = if codebook.min_dist == 4 {
            xi_miscall_4.iter().map(|xi| {
                LogProb::ln_sum_exp(&[
                    xi[(2, 1)], xi[(1, 2)], xi[(2, 3)], xi[(3, 2)]
                ])
            }).collect_vec()
        } else {
            vec![]
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
        event_probs.extend(prob_miscall_exact.iter().map(|&p| *Prob::from(p)));
        event_probs.extend(prob_miscall_mismatch.iter().map(|&p| *Prob::from(p)));

        let feature_type = if feature_record.expressed() {
            FeatureType::Expressed
        } else {
            FeatureType::NotExpressed
        };

        FeatureModel {
            feature_id: feature,
            event_probs: event_probs,
            neighbors: neighbors,
            min_dist: codebook.min_dist,
            counts: counts,
            event_counts: RefCell::new(Vec::new()),
            feature_type: feature_type
        }
    }

    /// Calculate maximum likelihood estimate of miscalls given expression and store in the given
    /// arrays.
    /// Returns the total absolute changes to estimated miscall counts.
    pub fn mle_miscalls(
        &self,
        expressions: &Expressions,
        miscalls_exact: &mut Miscalls,
        miscalls_mismatch: &mut Miscalls,
        rng: &mut rand::StdRng,
        debug: bool
    ) -> u32 {
        if !self.feature_type.is_expressed() {
            return 0;
        }

        let x = expressions[self.feature_id as usize];

        let mut total_change = 0;

        // We need to ensure that the sum of the estimated miscalls does not exceed x.
        // Hence, we iterate randomly over the neighbors such that we can simply stop once that
        // happens without causing a bias.
        let mut idx = (0..self.neighbors.len()).collect_vec();
        rng.shuffle(&mut idx);

        // by stochastic rounding, we ensure that on average \sum x*p_i = x
        // numeric rounding would not work if many of the probabilities are small enough to yield
        // expectations less than 1
        let mut stochastic_round = |v: f64| v.floor() as u32 + (rng.next_f64() <= v % 1.0) as u32;

        let mut rest = x;
        for i in idx {
            let n = self.neighbors[i];
            let miscall_exact = self.prob_miscall_exact(i) * x as f64;
            total_change += miscalls_exact.set(
                self.feature_id, n, stochastic_round(miscall_exact), &mut rest
            );

            if self.min_dist == 4 {
                let miscall_mismatch = self.prob_miscall_mismatch(i) * x as f64;
                let miscall_mismatch = stochastic_round(miscall_mismatch);
                let change = miscalls_mismatch.set(
                    self.feature_id, n, miscall_mismatch, &mut rest
                );
                total_change += change;
            }
        }
        if debug {
            debug!("total miscall mismatch={}", miscalls_mismatch.total_from(self.feature_id));
        }

        total_change
    }

    /// Calculate maximum likelihood estimate of expression given miscalls.
    /// Returns the absolute amount of change to the expression.
    pub fn mle_expression(
        &self,
        expressions: &mut Expressions,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls,
        feature_models: &[FeatureModel]
    ) -> i32 {
        let x_ = match self.feature_type {
            FeatureType::Expressed | FeatureType::Noise => {
                let calls_exact = self.counts.exact.saturating_sub(
                    miscalls_exact.total_to(self.feature_id)
                );
                let calls_mismatch = self.counts.mismatch.saturating_sub(
                    miscalls_mismatch.total_to(self.feature_id)
                );
                let event_count = calls_exact + calls_mismatch +
                                  miscalls_exact.total_from(self.feature_id) +
                                  miscalls_mismatch.total_from(self.feature_id);

                let x = event_count as f64 / (1.0 - self.event_probs[2]);

                cmp::max(
                    x.floor() as u32,
                    // or at least as many as we have events coming from this feature
                    // TODO does this allow x to decrease properly?
                    event_count
                )
            },
            FeatureType::NotExpressed => 0,
            // FeatureType::Noise => {
            //     // get remaining free miscalls of not expressed features
            //     // the best noise estimate should fill all of these
            //     // self.neighbors.iter().filter_map(|&n| {
            //     //     if let FeatureType::NotExpressed = feature_models[n].feature_type {
            //     //         let free_miscalls = |miscalls: &Miscalls| {
            //     //             miscalls.max_total_to(n) - miscalls.total_to(n) -
            //     //             miscalls.get(self.feature_id, n)
            //     //         };
            //     //
            //     //         free_miscalls(miscalls_exact) + free_miscalls(miscalls_mismatch)
            //     //     } else { None }
            //     // }
            //     let x: f64 = self.neighbors.iter().enumerate().filter_map(|(i, &n)| {
            //         if feature_models[n].feature_type.is_expressed() {
            //             let free_miscalls = |miscalls: &Miscalls| {
            //                 miscalls.max_total_to(n) - miscalls.total_to(n) -
            //                 miscalls.get(self.feature_id, n)
            //             };
            //             Some(
            //                 (
            //                     free_miscalls(miscalls_exact) as f64 +
            //                     free_miscalls(miscalls_mismatch) as f64
            //                 ) / (self.prob_miscall_exact(i) + self.prob_miscall_mismatch(i))
            //             )
            //         } else { None }
            //     }).sum::<f64>() / self.neighbors.len() as f64;
            //     // let x = free_miscalls as f64 / total_prob;
            //
            //     cmp::max(
            //         x.floor() as u32,
            //         miscalls_exact.total_from(self.feature_id) +
            //         miscalls_mismatch.total_from(self.feature_id)
            //     )
            // }
        };

        let x = expressions.get_mut(self.feature_id as usize).unwrap();

        let change = *x as i32 - x_ as i32;

        *x = x_;

        change
    }

    pub fn prob_dropout(&self) -> f64 {
        self.event_probs[2]
    }

    pub fn prob_miscall_exact(&self, neighbor_index: usize) -> f64 {
        self.event_probs[3 + neighbor_index]
    }

    pub fn prob_miscall_mismatch(&self, neighbor_index: usize) -> f64 {
        self.event_probs[3 + self.neighbors.len() + neighbor_index]
    }

    /// Minimum expression given fixed miscalls.
    pub fn min_expression(&self, miscalls_exact: &Miscalls, miscalls_mismatch: &Miscalls) -> u32 {
        self.counts.exact.saturating_sub(miscalls_exact.total_to(self.feature_id)) +
        self.counts.mismatch.saturating_sub(miscalls_mismatch.total_to(self.feature_id)) +
        miscalls_exact.total_from(self.feature_id) + miscalls_mismatch.total_from(self.feature_id)
    }

    /// Calculate likelihood of given expression.
    pub fn likelihood(
        &self,
        x: u32,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> LogProb {
        let count_exact_miscalled = miscalls_exact.total_to(self.feature_id);
        let count_mismatch_miscalled = miscalls_mismatch.total_to(self.feature_id);

        if self.counts.exact < count_exact_miscalled {
            panic!("bug: likelihood cannot be calculated if exact miscalls exceed counts: counts={}, miscalls={}", self.counts.exact, count_exact_miscalled);
        }
        if self.counts.mismatch < count_mismatch_miscalled {
            panic!("bug: likelihood cannot be calculated if mismatch miscalls exceed counts: counts={}, miscalls={}", self.counts.mismatch, count_mismatch_miscalled);
        }
        let total_count = self.counts.exact - count_exact_miscalled +
                          self.counts.mismatch - count_mismatch_miscalled;
        if total_count > x {
            panic!("bug: event count exceeds expression: x={}, event count={}", x, total_count);
        }

        let calls_exact = self.counts.exact - count_exact_miscalled;
        let calls_mismatch = self.counts.mismatch - count_mismatch_miscalled;

        // setup event counts
        {
            let mut event_counts = self.event_counts.borrow_mut();
            event_counts.clear();
            event_counts.push(
                // exact calls
                calls_exact
            );
            event_counts.push(
                // mismatch calls
                calls_mismatch
            );
            event_counts.push(
                // dropouts (fill in later)
                0
            );
            // miscalls
            event_counts.extend(self.neighbors.iter().map(
                |&n| miscalls_exact.get(self.feature_id, n)
            ));
            if self.min_dist == 4 {
                event_counts.extend(
                    self.neighbors.iter().map(
                        |&n| miscalls_mismatch.get(self.feature_id, n)
                    )
                );
            }

            // set dropouts
            let total_counts = event_counts.iter().sum();
            assert!(x >= total_counts, "total event counts exceed expression: x={} {:?}", x, event_counts);
            let dropouts = x - total_counts;
            event_counts[2] = dropouts;
        }

        LogProb::from(Prob(
            multinomial_pdf(&self.event_probs, &self.event_counts.borrow())
        ))
    }
}
