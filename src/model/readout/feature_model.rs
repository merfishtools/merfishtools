use std::cell::RefCell;
use std::cmp;

use rand;
use rand::Rng;
use itertools::Itertools;
use rgsl::randist::multinomial::multinomial_pdf;

use bio::stats::{Prob, LogProb};

use io::codebook::{Codebook, FeatureID};

use model::readout::{Expressions, Counts, Miscalls, Xi};


pub trait AbstractFeatureModel {
    /// FeatureID of this feature.
    fn feature_id(&self) -> FeatureID;

    /// Neighbors of this feature.
    fn neighbors(&self) -> &[FeatureID];

    /// Minimum hamming distance between any pair of features.
    fn min_dist(&self) -> u8;

    /// Number of events (calls and miscalls) except dropout.
    fn event_count_non_dropout(
        &self,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> u32;

    /// Dropout probability.
    fn prob_dropout(&self) -> f64;

    /// Probability of exact miscall to given neighbor index.
    fn prob_miscall_exact(&self, neighbor_index: usize) -> f64;

    /// Probability of mismatching miscall to given neighbor index.
    fn prob_miscall_mismatch(&self, neighbor_index: usize) -> f64;

    /// Calculate maximum likelihood estimate of miscalls given expression and store in the given
    /// arrays.
    /// Returns the total absolute changes to estimated miscall counts.
    fn mle_miscalls(
        &self,
        expressions: &Expressions,
        miscalls_exact: &mut Miscalls,
        miscalls_mismatch: &mut Miscalls,
        rng: &mut rand::StdRng
    ) {
        let x = expressions[self.feature_id()];

        // We need to ensure that the sum of the estimated miscalls does not exceed x.
        // Hence, we iterate randomly over the neighbors such that we can simply stop once that
        // happens without causing a bias.
        let mut idx = (0..self.neighbors().len()).collect_vec();
        rng.shuffle(&mut idx);

        // by stochastic rounding, we ensure that on average \sum x*p_i = x
        // numeric rounding would not work if many of the probabilities are small enough to yield
        // expectations less than 1
        let mut stochastic_round = |v: f64| v.floor() as u32 + (rng.next_f64() <= v % 1.0) as u32;

        let mut rest = x;
        for i in idx {
            let n = self.neighbors()[i];
            let miscall_exact = self.prob_miscall_exact(i) * x as f64;
            miscalls_exact.set(
                self.feature_id(), n, stochastic_round(miscall_exact), &mut rest
            );

            if self.min_dist() == 4 {
                let miscall_mismatch = self.prob_miscall_mismatch(i) * x as f64;
                let miscall_mismatch = stochastic_round(miscall_mismatch);
                miscalls_mismatch.set(
                    self.feature_id(), n, miscall_mismatch, &mut rest
                );
            }
        }
    }

    /// Calculate maximum likelihood estimate of expression given miscalls.
    /// Returns the absolute amount of change to the expression.
    fn mle_expression(
        &self,
        expressions: &mut Expressions,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> i32 {
        let event_count = self.event_count_non_dropout(miscalls_exact, miscalls_mismatch);

        let x_ = cmp::max(
            (event_count as f64 / (1.0 - self.prob_dropout())).floor() as u32,
            // or at least as many as we have events coming from this feature
            // TODO does this allow x to decrease properly?
            event_count
        );

        let x = expressions.get_mut(self.feature_id()).unwrap();

        let change = *x as i32 - x_ as i32;

        *x = x_;

        change
    }
}


pub struct FeatureModel {
    pub(crate) feature_id: FeatureID,
    event_probs: Vec<f64>,
    neighbors: Vec<FeatureID>,
    min_dist: u8,
    pub(crate) counts: Counts,
    event_counts: RefCell<Vec<u32>>
}


impl FeatureModel {
    pub fn new(
        feature: FeatureID,
        counts: Counts,
        codebook: &Codebook,
        xi: &Xi
    ) -> Self {
        assert!(
            codebook.min_dist == 2 || codebook.min_dist == 4,
            "unsupported hamming distance: {}, supported: 2, 4",
            codebook.min_dist
        );

        let feature_record = codebook.record(feature);

        let xi_call = xi.prob(feature_record.codeword(), feature_record.codeword());

        let prob_call_exact = xi_call[0];
        let prob_call_mismatch = if codebook.min_dist == 4 {
            xi_call[1]
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
            xi.prob(feature_record.codeword(), codebook.record(n).codeword())
        }).collect_vec();
        let xi_miscall_2 = if codebook.min_dist == 2 {
            codebook.neighbors(feature, 2).iter().map(|&n| {
                xi.prob(feature_record.codeword(), codebook.record(n).codeword())
            }).collect_vec()
        } else {
            vec![]
        };

        let prob_miscall_exact = xi_miscall_4.iter().map(|xi| {
            xi[0]
        }).chain(xi_miscall_2.iter().map(|xi| xi[0])).collect_vec();

        let prob_miscall_mismatch = if codebook.min_dist == 4 {
            xi_miscall_4.iter().map(|xi| xi[1]).collect_vec()
        } else {
            vec![]
        };

        let mut total = vec![prob_call_exact, prob_call_mismatch];
        total.extend(&prob_miscall_exact);
        total.extend(&prob_miscall_mismatch);

        let prob_total = LogProb::ln_sum_exp(&total);

        let mut event_probs = vec![
            *Prob::from(prob_call_exact),
            *Prob::from(prob_call_mismatch),
            *Prob::from(prob_total.ln_one_minus_exp())
        ];
        event_probs.extend(prob_miscall_exact.iter().map(|&p| *Prob::from(p)));
        event_probs.extend(prob_miscall_mismatch.iter().map(|&p| *Prob::from(p)));

        debug!("{:?}", event_probs);

        FeatureModel {
            feature_id: feature,
            event_probs: event_probs,
            neighbors: neighbors,
            min_dist: codebook.min_dist,
            counts: counts,
            event_counts: RefCell::new(Vec::new())
        }
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

    /// Minimum expression given fixed miscalls.
    pub fn min_expression(&self, miscalls_exact: &Miscalls, miscalls_mismatch: &Miscalls) -> u32 {
        self.counts.exact.saturating_sub(miscalls_exact.total_to(self.feature_id)) +
        self.counts.mismatch.saturating_sub(miscalls_mismatch.total_to(self.feature_id)) +
        miscalls_exact.total_from(self.feature_id) + miscalls_mismatch.total_from(self.feature_id)
    }
}


impl AbstractFeatureModel for FeatureModel {
    fn feature_id(&self) -> FeatureID {
        self.feature_id
    }

    fn neighbors(&self) -> &[FeatureID] {
        &self.neighbors
    }

    fn min_dist(&self) -> u8 {
        self.min_dist
    }

    fn event_count_non_dropout(
        &self,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> u32 {
        let calls_exact = self.counts.exact.saturating_sub(
            miscalls_exact.total_to(self.feature_id)
        );
        let calls_mismatch = self.counts.mismatch.saturating_sub(
            miscalls_mismatch.total_to(self.feature_id)
        );

        calls_exact + calls_mismatch +
        miscalls_exact.total_from(self.feature_id) +
        miscalls_mismatch.total_from(self.feature_id)
    }

    fn prob_dropout(&self) -> f64 {
        self.event_probs[2]
    }

    fn prob_miscall_exact(&self, neighbor_index: usize) -> f64 {
        self.event_probs[3 + neighbor_index]
    }

    fn prob_miscall_mismatch(&self, neighbor_index: usize) -> f64 {
        self.event_probs[3 + self.neighbors.len() + neighbor_index]
    }
}



pub struct NoiseModel {
    pub(crate) feature_id: FeatureID,
    pub(crate) not_expressed_feature_ids: Vec<FeatureID>,
    pub(crate) not_expressed_counts: Vec<Counts>,
    probs_miscall_exact: Vec<Prob>,
    probs_miscall_mismatch: Vec<Prob>,
    prob_dropout: Prob,
    neighbors: Vec<FeatureID>,
    min_dist: u8,
}


impl NoiseModel {
    pub fn new(
        feature_id: FeatureID,
        not_expressed_feature_ids: Vec<FeatureID>,
        not_expressed_counts: Vec<Counts>,
        codebook: &Codebook,
        xi: &Xi
    ) -> Self {
        assert!(
            codebook.min_dist == 2 || codebook.min_dist == 4,
            "unsupported hamming distance: {}, supported: 2, 4",
            codebook.min_dist
        );

        // miscalls generated from noise
        let xi_miscall = codebook.records().iter().filter_map(|feat_rec| {
            if feat_rec.expressed() {
                Some(xi.prob(codebook.noise().codeword(), feat_rec.codeword()))
            } else { None }
        }).collect_vec();
        let neighbors = codebook.records().iter().filter_map(|rec| {
            if rec.expressed() {
                Some(codebook.get_id(rec.name()))
            } else { None }
        }).collect_vec();

        let probs_miscall_exact = xi_miscall.iter().map(|xi| xi[0]).collect_vec();
        let probs_miscall_mismatch = if codebook.min_dist == 4 {
            xi_miscall.iter().map(|xi| xi[1]).collect_vec()
        } else {
            vec![]
        };

        let prob_total_miscall_exact = LogProb::ln_sum_exp(&probs_miscall_exact);
        let prob_total_miscall_mismatch = LogProb::ln_sum_exp(&probs_miscall_mismatch);

        // probs for calls of not expressed features
        let xi_not_expressed = not_expressed_feature_ids.iter().map(|&f| {
            xi.prob(codebook.noise().codeword(), codebook.record(f).codeword())
        }).collect_vec();

        let prob_not_expressed_exact = LogProb::ln_sum_exp(
            &xi_not_expressed.iter().map(|xi| xi[0]).collect_vec()
        );
        let prob_not_expressed_mismatch = LogProb::ln_sum_exp(
            &xi_not_expressed.iter().map(|xi| xi[1]).collect_vec()
        );

        // dropout
        let prob_dropout = LogProb::ln_sum_exp(&[
            prob_not_expressed_exact,
            prob_not_expressed_mismatch,
            prob_total_miscall_exact,
            prob_total_miscall_mismatch
        ]).ln_one_minus_exp();

        NoiseModel {
            feature_id: feature_id,
            not_expressed_feature_ids: not_expressed_feature_ids,
            not_expressed_counts: not_expressed_counts,
            probs_miscall_exact: probs_miscall_exact.iter().map(|&p| Prob::from(p)).collect_vec(),
            probs_miscall_mismatch: probs_miscall_mismatch.iter().map(|&p| Prob::from(p)).collect_vec(),
            prob_dropout: Prob::from(prob_dropout),
            neighbors: neighbors,
            min_dist: codebook.min_dist
        }
    }
}


impl AbstractFeatureModel for NoiseModel {
    fn feature_id(&self) -> FeatureID {
        self.feature_id
    }

    fn neighbors(&self) -> &[FeatureID] {
        &self.neighbors
    }

    fn min_dist(&self) -> u8 {
        self.min_dist
    }

    fn event_count_non_dropout(
        &self,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> u32 {
        let mut calls_exact = 0;
        let mut calls_mismatch = 0;
        for (&id, counts) in self.not_expressed_feature_ids.iter().zip(&self.not_expressed_counts) {
            calls_exact += counts.exact.saturating_sub(miscalls_exact.total_to(id));
            calls_mismatch += counts.mismatch.saturating_sub(miscalls_mismatch.total_to(id));
        }

        let event_count = calls_exact + calls_mismatch +
                          miscalls_exact.total_from(self.feature_id) +
                          miscalls_mismatch.total_from(self.feature_id);

        event_count
    }

    fn prob_dropout(&self) -> f64 {
        *self.prob_dropout
    }

    fn prob_miscall_exact(&self, neighbor_index: usize) -> f64 {
        *self.probs_miscall_exact[neighbor_index]
    }

    fn prob_miscall_mismatch(&self, neighbor_index: usize) -> f64 {
        *self.probs_miscall_mismatch[neighbor_index]
    }
}
