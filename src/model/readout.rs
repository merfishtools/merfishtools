use std::mem;
use std::cell::RefCell;

use itertools::Itertools;
use ndarray::prelude::*;
use statrs::distribution::{Discrete, Multinomial};
use bit_vec::BitVec;

use bio::stats::{Prob, LogProb};

use io::codebook::{self, Codebook, Codeword, FeatureID};


const MAX_ERR: u8 = 6;


pub type Expressions = Vec<u32>;


pub struct JointModel {
    feature_models: Vec<FeatureModel>,
    expressions: Expressions,
    miscalls: Miscalls,
    margin: u32
}


impl JointModel {
    pub fn new<'a, I: Iterator<Item=(&'a String, &'a Counts)>>(
        counts: I,
        p0: &[Prob],
        p1: &[Prob],
        codebook: &Codebook,
        window_width: u32) -> Self {
        let mut xi = Xi::new(p0, p1, MAX_ERR);

        let feature_models = counts.map(
            |(feature, counts)| {
                FeatureModel::new(feature, counts.clone(), codebook, &xi)
            }
        ).collect_vec();

        // calculate start values (we take raw counts as start expression)
        let expressions = feature_models.iter().map(
            |model| model.counts.total()
        ).collect_vec();

        let miscalls = Miscalls::new(feature_models.len());

        JointModel {
            feature_models: feature_models,
            expressions: expressions,
            miscalls: miscalls,
            margin: window_width / 2
        }
    }

    /// Jointly estimate MAP using the EM algorithm.
    ///
    /// Expressions are our model parameters to estimate, miscalls are our latent variables,
    /// counts are our observed variables.
    pub fn expectation_maximization(&mut self) {
        // EM iterations
        loop {
            let total: u32 = self.expressions.iter().sum::<u32>() + self.miscalls.sum();

            let mut total_change: u32 = 0;
            // E-step: estimate miscalls
            for m in &self.feature_models {
                total_change += m.mle_miscalls(&self.expressions, &mut self.miscalls);
            }

            // M-step: estimate expressions that maximize the probability
            for m in &self.feature_models {
                total_change += m.mle_expression(&mut self.expressions, &self.miscalls)
            }

            // check for convergence (less than 1% change)
            // TODO tune, parameterize, max number of iterations?
            if total_change as f64 / total as f64 <= 0.01 {
                break;
            }
        }
    }

    pub fn likelihood(&mut self, feature_id: FeatureID, x: u32) -> LogProb {
        self.feature_models[feature_id as usize].likelihood(x, &self.miscalls)
    }

    pub fn map_estimate(&self, feature_id: FeatureID) -> u32 {
        self.expressions[feature_id as usize]
    }

    /// Prior window for calculating expression PMF
    pub fn window(&self, feature_id: FeatureID) -> (u32, u32) {
        let x = self.map_estimate(feature_id);

        (self.feature_models[feature_id as usize].min_expression(&self.miscalls), x + self.margin)
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


pub struct Miscalls {
    exact: Array2<u32>,
    mismatch: Array2<u32>
}


impl Miscalls {
    /// Create new instance.
    pub fn new(feature_count: usize) -> Self {
        Miscalls {
            exact: Array::from_elem((feature_count, feature_count), 0),
            mismatch: Array::from_elem((feature_count, feature_count), 0)
        }
    }

    pub fn exact(&self, i: FeatureID, j: FeatureID) -> u32 {
        self.exact[(i as usize, j as usize)]
    }

    pub fn mismatch(&self, i: FeatureID, j: FeatureID) -> u32 {
        self.mismatch[(i as usize, j as usize)]
    }

    pub fn set_exact(&mut self, i: FeatureID, j: FeatureID, value: u32) -> u32 {
        let entry = self.exact.get_mut((i as usize, j as usize)).unwrap();
        let change = (*entry as i32 - value as i32).abs() as u32;
        *entry = value;

        change
    }

    pub fn set_mismatch(&mut self, i: FeatureID, j: FeatureID, value: u32) -> u32 {
        let entry = self.mismatch.get_mut((i as usize, j as usize)).unwrap();
        let change = (*entry as i32 - value as i32).abs() as u32;
        *entry = value;

        change
    }

    /// Total miscalls coming from given feature as the true origin.
    pub fn total_exact_from(&self, feature: FeatureID) -> u32 {
        self.exact.row(feature).scalar_sum()
    }

    /// Total miscalls coming from given feature as the true origin.
    pub fn total_mismatch_from(&self, feature: FeatureID) -> u32 {
        self.mismatch.row(feature).scalar_sum()
    }

    /// Total miscalls leading to given feature.
    pub fn total_exact_to(&self, feature: FeatureID) -> u32 {
        self.exact.column(feature).scalar_sum()
    }

    /// Total miscalls leading to given feature.
    pub fn total_mismatch_to(&self, feature: FeatureID) -> u32 {
        self.mismatch.column(feature).scalar_sum()
    }

    pub fn sum(&self) -> u32 {
        self.exact.scalar_sum() + self.mismatch.scalar_sum()
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
        let xi_miscall_2 = codebook.neighbors(feature, 2).iter().map(|&n| {
            xi.prob_error(feature_record, Some(codebook.get_record(n)))
        }).collect_vec();

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
    pub fn mle_miscalls(&self, expressions: &Expressions, miscalls: &mut Miscalls) -> u32 {
        let x = expressions[self.feature_id as usize];
        let expectation = |p: f64| (p * x as f64).round() as u32;

        let mut total_change = 0;

        let offset_exact = 3;
        let offset_mismatch = 3 + self.neighbors.len();
        for &n in &self.neighbors {
            let prob_miscall_exact = self.event_probs[offset_exact + n as usize];
            let prob_miscall_mismatch = self.event_probs[offset_mismatch + n as usize];

            total_change += miscalls.set_exact(
                self.feature_id, n, expectation(prob_miscall_exact)
            );
            total_change += miscalls.set_mismatch(
                self.feature_id, n, expectation(prob_miscall_mismatch)
            );
        }

        total_change
    }

    /// Calculate maximum likelihood estimate of expression given miscalls.
    /// Returns the absolute amount of change to the expression.
    pub fn mle_expression(&self, expressions: &mut Expressions, miscalls: &Miscalls) -> u32 {
        let x1 = (self.counts.exact - miscalls.total_exact_to(self.feature_id)) as f64 /
                 self.event_probs[0];
        let x2 = if self.min_dist == 4 {
            (self.counts.mismatch - miscalls.total_mismatch_to(self.feature_id)) as f64 /
             self.event_probs[1]
        } else {
            x1
        };

        let x_ = (x1 as f64 + x2 as f64 / 2.0).round() as u32;
        let x = expressions.get_mut(self.feature_id as usize).unwrap();

        let change = (*x as i32 - x_ as i32).abs() as u32;

        *x = x_;

        change
    }

    /// Minimum expression given fixed miscalls.
    pub fn min_expression(&self, miscalls: &Miscalls) -> u32 {
        self.counts.exact.saturating_sub(miscalls.total_exact_to(self.feature_id)) +
        self.counts.mismatch.saturating_sub(miscalls.total_mismatch_to(self.feature_id)) +
        miscalls.total_exact_from(self.feature_id) + miscalls.total_mismatch_from(self.feature_id)
    }

    /// Calculate likelihood of given expression.
    pub fn likelihood(
        &self,
        x: u32,
        miscalls: &Miscalls
    ) -> LogProb {
        let count_exact_miscalled = miscalls.total_exact_to(self.feature_id);
        let count_mismatch_miscalled = miscalls.total_mismatch_to(self.feature_id);

        if self.counts.exact < count_exact_miscalled {
            return LogProb::ln_zero();
        }
        if self.counts.mismatch < count_mismatch_miscalled {
            return LogProb::ln_zero();
        }
        if (self.counts.exact - count_exact_miscalled + self.counts.mismatch - count_mismatch_miscalled) > x {
            return LogProb::ln_zero();
        }

        // setup event counts
        {
            let mut event_counts = self.event_counts.borrow_mut();
            event_counts.clear();
            event_counts.push(
                // exact calls
                (self.counts.exact - count_exact_miscalled) as u64
            );
            event_counts.push(
                // mismatch calls
                (self.counts.mismatch - count_mismatch_miscalled) as u64
            );
            event_counts.push(
                // dropouts
                x.saturating_sub(self.counts.total()) as u64
            );
            // miscalls
            event_counts.extend(self.neighbors.iter().map(
                |&n| miscalls.exact(self.feature_id, n) as u64
            ));
            if self.min_dist == 4 {
                event_counts.extend(
                    self.neighbors.iter().map(|&n| miscalls.mismatch(self.feature_id, n) as u64)
                );
            }
        }

        // check if sum matches x, else return 0
        if self.event_counts.borrow().iter().sum::<u64>() != x as u64 {
            debug!("Likelihood calculated for unrealistic expression. This should be optimized.");
            return LogProb::ln_zero();
        }

        let multinomial = Multinomial::new(&self.event_probs, x as u64).unwrap();

        LogProb::from(Prob(multinomial.pmf(&self.event_counts.borrow())))
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
