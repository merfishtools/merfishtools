use std::mem;
use std::cell::RefCell;
use std::cmp;

use rand;
use rand::Rng;
use itertools::Itertools;
use ndarray::prelude::*;
use statrs::distribution::{Discrete, Multinomial, Distribution};
use bit_vec::BitVec;

use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codebook, Codeword, FeatureID};


const MAX_ERR: u8 = 6;


pub type Expressions = Vec<u32>;


pub struct JointModel {
    feature_models: Vec<FeatureModel>,
    expressions: Expressions,
    miscalls_exact: Miscalls,
    miscalls_mismatch: Miscalls,
    margin: u32,
    em_run: bool,
    rng: rand::StdRng
}


impl JointModel {
    pub fn new<'a, I: Iterator<Item=(&'a String, &'a Counts)>>(
        counts: I,
        p0: &[Prob],
        p1: &[Prob],
        codebook: &Codebook,
        window_width: u32) -> Self {
        let mut xi = Xi::new(p0, p1, MAX_ERR);

        let mut feature_models = counts.map(
            |(feature, counts)| {
                FeatureModel::new(feature, counts.clone(), codebook, &xi)
            }
        ).collect_vec();
        // sort feature models by id, such that we can access them by id in the array
        feature_models.sort_by_key(|m| m.feature_id);

        // calculate start values (we take raw counts as start expression)
        let expressions = feature_models.iter().map(
            |model| model.counts.total()
        ).collect_vec();

        let miscalls_exact = Miscalls::new(&feature_models, true);
        let miscalls_mismatch = Miscalls::new(&feature_models, false);

        for (i, m) in feature_models.iter().enumerate() {
            assert_eq!(i, m.feature_id);
        }

        JointModel {
            feature_models: feature_models,
            expressions: expressions,
            miscalls_exact: miscalls_exact,
            miscalls_mismatch: miscalls_mismatch,
            margin: window_width / 2,
            em_run: false,
            rng: rand::StdRng::new().unwrap()
        }
    }

    /// Jointly estimate MAP using the EM algorithm.
    ///
    /// Expressions are our model parameters to estimate, miscalls are our latent variables,
    /// counts are our observed variables.
    pub fn expectation_maximization(&mut self) {
        let mut idx = (0..self.feature_models.len()).collect_vec();
        let mut phase = 0;
        let mut count_no_change = 0;
        let mut map_likelihood = LogProb::ln_zero();
        let mut map_miscalls_exact = None;
        let mut map_miscalls_mismatch = None;
        let mut map_expressions = None;
        let mut i = 0;
        // EM iterations
        loop {
            // Shuffle models, such that miscalls can be drawn unbiased
            // this is necessary because there is a maximum miscall count for each feature.
            // If we would not shuffle, it would be easier for the first model to provide miscalls.
            self.rng.shuffle(&mut idx);

            let total: u32 = self.expressions.iter().sum::<u32>();

            let mut e_change: u32 = 0;
            let mut m_change: u32 = 0;
            // E-step: estimate miscalls
            for &j in &idx {
                e_change += self.feature_models[j].mle_miscalls(
                    &self.expressions,
                    &mut self.miscalls_exact,
                    &mut self.miscalls_mismatch,
                    &mut self.rng
                );
            }

            // M-step: estimate expressions that maximize the probability
            for m in &self.feature_models {
                m_change += m.mle_expression(
                    &mut self.expressions,
                    &self.miscalls_exact,
                    &self.miscalls_mismatch
                );
            }

            let fraction_change = m_change as f64 / total as f64;

            if phase == 0 {
                // phase 0: run until convergence
                if fraction_change <= 0.05 {
                    count_no_change += 1;
                } else {
                    count_no_change = 0;
                }
                // convergence or at least 100 steps
                if count_no_change >= 5 || i >= 100 {
                    phase = 1;
                }
                debug!("EM-iteration (phase 0): e-change={}, m-change={}, %={}",
                    e_change,
                    m_change,
                    fraction_change,
                );

            } else if phase == 1 {
                // phase 1: take best likelihood from 100 iterations
                let mut likelihood = LogProb(0.0);
                for m in &self.feature_models {
                    let x = self.expressions[m.feature_id];
                    let l = m.likelihood(x, &self.miscalls_exact, &self.miscalls_mismatch);
                    likelihood = likelihood + l;
                }

                if likelihood > map_likelihood {
                    map_miscalls_exact = Some(self.miscalls_exact.clone());
                    map_miscalls_mismatch = Some(self.miscalls_mismatch.clone());
                    map_expressions = Some(self.expressions.clone());
                }

                debug!(
                    "EM-iteration (phase 1): e-change={}, m-change={}, %={}, L={}",
                    e_change,
                    m_change,
                    fraction_change,
                    *likelihood
                );

                if i >= 100 {
                    break;
                }
            }

            i += 1;
        }
        self.miscalls_exact = map_miscalls_exact.unwrap();
        self.miscalls_mismatch = map_miscalls_mismatch.unwrap();
        self.expressions = map_expressions.unwrap();
        self.em_run = true;
    }

    pub fn likelihood(&mut self, feature_id: FeatureID, x: u32) -> LogProb {
        assert!(self.em_run);
        self.feature_models[feature_id as usize].likelihood(
            x, &self.miscalls_exact, &self.miscalls_mismatch
        )
    }

    pub fn map_estimate(&self, feature_id: FeatureID) -> u32 {
        assert!(self.em_run);
        self.expressions[feature_id as usize]
    }

    pub fn feature_model(&self, feature_id: FeatureID) -> &FeatureModel {
        &self.feature_models[feature_id]
    }

    /// Prior window for calculating expression PMF
    pub fn window(&self, feature_id: FeatureID) -> (u32, u32) {
        assert!(self.em_run);
        let x = self.map_estimate(feature_id);
        let xmin = self.feature_models[feature_id as usize].min_expression(
            &self.miscalls_exact, &self.miscalls_mismatch
        );
        let xmax = x + self.margin;
        assert!(xmax >= xmin, "xmax has to be greater than xmin: xmin={}, xmax={}", xmin, xmax);
        (xmin, xmax)
    }
}


#[derive(Debug, Clone)]
pub struct Counts {
    pub exact: u32,
    pub mismatch: u32
}


impl Counts {
    pub fn total(&self) -> u32 {
        self.exact + self.mismatch
    }
}


#[derive(Clone)]
pub struct Miscalls {
    miscalls: Array2<u32>,
    total_to: Vec<u32>,
    max_total_to: Vec<u32>,
    total: u32
}


impl Miscalls {
    /// Create new instance.
    pub fn new(feature_models: &[FeatureModel], exact: bool) -> Self {
        let feature_count = feature_models.len();
        let mut max_total_to = vec![0; feature_models.len()];
        for m in feature_models.iter() {
            max_total_to[m.feature_id] = if exact { m.counts.exact } else { m.counts.mismatch };
        }
        Miscalls {
            miscalls: Array::from_elem((feature_count, feature_count), 0),
            total_to: vec![0; feature_count],
            max_total_to: max_total_to,
            total: 0
        }
    }

    pub fn get(&self, i: FeatureID, j: FeatureID) -> u32 {
        self.miscalls[(i as usize, j as usize)]
    }

    pub fn set(&mut self, i: FeatureID, j: FeatureID, value: u32) -> u32 {
        let change = {
            let entry = self.miscalls.get_mut((i as usize, j as usize)).unwrap();
            let mut change = value as i32 - *entry as i32;

            // update column sum
            let colsum = self.total_to[j] as i32 + change;
            // get difference to maximum, only consider exceeding
            let maxdiff = cmp::max(colsum - self.max_total_to[j] as i32, 0);
            // correct change
            change -= maxdiff;
            self.total_to[j] = (colsum - maxdiff) as u32;

            // update total
            self.total = (self.total as i32 + change) as u32;
            // update entry
            *entry = (*entry as i32 + change) as u32;

            change.abs() as u32
        };
        assert!(self.total_to[j] <= self.max_total_to[j]);
        assert_eq!(self.total_to[j], self.miscalls.column(j).scalar_sum());

        change
    }

    /// Total miscalls coming from given feature as the true origin.
    pub fn total_from(&self, feature: FeatureID) -> u32 {
        self.miscalls.row(feature).scalar_sum()
    }

    /// Total miscalls leading to given feature.
    pub fn total_to(&self, feature: FeatureID) -> u32 {
        self.total_to[feature]
    }

    pub fn total(&self) -> u32 {
        self.total
    }
}


pub struct FeatureModel {
    feature_id: FeatureID,
    event_probs: Vec<f64>,
    neighbors: Vec<FeatureID>,
    min_dist: u8,
    counts: Counts,
    event_counts: RefCell<Vec<u64>>
}


impl FeatureModel {
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

        let xi_miscall_4 = codebook.neighbors(feature, 4).iter().map(|&n| {
            xi.prob_error(feature_record, Some(codebook.get_record(n)))
        }).collect_vec();
        let xi_miscall_2 = if codebook.min_dist == 2 {
            codebook.neighbors(feature, 2).iter().map(|&n| {
                xi.prob_error(feature_record, Some(codebook.get_record(n)))
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

        FeatureModel {
            feature_id: feature,
            event_probs: event_probs,
            neighbors: neighbors,
            min_dist: codebook.min_dist,
            counts: counts,
            event_counts: RefCell::new(Vec::new())
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
        rng: &mut rand::StdRng
    ) -> u32 {
        let x = expressions[self.feature_id as usize];

        let multinomial = Multinomial::new(&self.event_probs, x as u64).unwrap();
        let sample = multinomial.sample(rng);

        let offset_exact = 3;
        let offset_mismatch = 3 + self.neighbors.len();
        let mut total_change = 0;

        for (i, &n) in self.neighbors.iter().enumerate() {
            let miscall_exact = sample[offset_exact + i] as u32;
            let miscall_mismatch = sample[offset_mismatch + i] as u32;

            total_change += miscalls_exact.set(
                self.feature_id, n, miscall_exact
            );
            total_change += miscalls_mismatch.set(
                self.feature_id, n, miscall_mismatch
            );
        }

        total_change
    }

    /// Calculate maximum likelihood estimate of expression given miscalls.
    /// Returns the absolute amount of change to the expression.
    pub fn mle_expression(
        &self,
        expressions: &mut Expressions,
        miscalls_exact: &Miscalls,
        miscalls_mismatch: &Miscalls
    ) -> u32 {
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

        let x_ = cmp::max(
            x.floor() as u32,
            // or at least as many as we have events coming from this feature
            event_count
        );

        let x = expressions.get_mut(self.feature_id as usize).unwrap();

        let change = (*x as i32 - x_ as i32).abs() as u32;

        *x = x_;

        change
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
                calls_exact as u64
            );
            event_counts.push(
                // mismatch calls
                calls_mismatch as u64
            );
            event_counts.push(
                // dropouts (fill in later)
                0
            );
            // miscalls
            event_counts.extend(self.neighbors.iter().map(
                |&n| miscalls_exact.get(self.feature_id, n) as u64
            ));
            if self.min_dist == 4 {
                event_counts.extend(
                    self.neighbors.iter().map(
                        |&n| miscalls_mismatch.get(self.feature_id, n) as u64
                    )
                );
            }

            // set dropouts
            let total_counts = event_counts.iter().sum();
            assert!(x as u64 >= total_counts, "total event counts exceed expression: x={} {:?}", x, event_counts);
            let dropouts = x as u64 - total_counts;
            event_counts[2] = dropouts;
        }

        let multinomial = Multinomial::new(&self.event_probs, x as u64).unwrap();

        LogProb(multinomial.ln_pmf(&self.event_counts.borrow()))
    }
}


/// Basic model for the probability of making i 1-0 and j 0-1 errors.
pub struct Xi {
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

    fn prob_error(
        &self,
        source: &codebook::Record,
        target: Option<&codebook::Record>
    ) -> Array2<LogProb> {
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
