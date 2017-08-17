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

use model::readout::{Expressions, Counts, Miscalls, FeatureModel, Xi};


pub struct JointModel {
    feature_models: Vec<FeatureModel>,
    expressions: Expressions,
    miscalls_exact: Miscalls,
    miscalls_mismatch: Miscalls,
    margin: u32,
    em_run: bool,
    rng: rand::StdRng,
    debug_id: Option<FeatureID>,
    noise_id: FeatureID
}


impl JointModel {
    pub fn new<'a, I: Iterator<Item=(&'a String, &'a Counts)>>(
        counts: I,
        p0: &[Prob],
        p1: &[Prob],
        codebook: &Codebook,
        window_width: u32) -> Self {
        let mut xi = Xi::new(p0, p1);
        let mut rng = rand::StdRng::new().unwrap();

        let mut feature_models = counts.map(
            |(feature, counts)| {
                FeatureModel::new(feature, counts.clone(), codebook, &xi)
            }
        ).collect_vec();
        let noise_id = feature_models.len();
        feature_models.push(FeatureModel::noise(noise_id, codebook, &xi));
        // sort feature models by id, such that we can access them by id in the array
        feature_models.sort_by_key(|m| m.feature_id);


        // calculate start values (we take raw counts as start expression)
        let mut expressions = Array1::from_iter(feature_models.iter().map(
            |model| if model.feature_type.is_expressed() { rng.next_u32() } else { eprintln!("misid={}", model.counts.total()); 0 }
        ));
        // noise start value is sum of total expressions
        expressions[noise_id] = expressions.scalar_sum();

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
            rng: rng,
            debug_id: None, //codebook.get_id("THBS1")
            noise_id: noise_id
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
        // let mut map_likelihood = LogProb::ln_zero();
        // let mut map_miscalls_exact = None;
        // let mut map_miscalls_mismatch = None;
        // let mut map_expressions = None;
        let mut last_expressions = self.expressions.clone();
        let n_iterations = 20;
        let mut i = 0;

        let mut last_changes = Array2::from_elem((self.expressions.len(), 5), 0.0);

        // EM iterations
        loop {
            let total: u32 = self.expressions.len() as u32;
            //let total: u32 = self.expressions.iter().sum::<u32>();

            // E-step: estimate miscalls
            debug!("E-step");

            // Shuffle models, such that miscalls can be drawn unbiased
            // this is necessary because there is a maximum miscall count for each feature.
            // If we would not shuffle, it would be easier for the first model to provide miscalls.
            self.rng.shuffle(&mut idx);
            for &j in &idx {
                self.feature_models[j].mle_miscalls(
                    &self.expressions,
                    &mut self.miscalls_exact,
                    &mut self.miscalls_mismatch,
                    &mut self.rng,
                    j == self.noise_id
                );
            }

            if i == n_iterations {
                break;
            }

            // M-step: estimate expressions that maximize the probability
            debug!("M-step");
            for (i, m) in self.feature_models.iter().enumerate() {
                last_changes[(i, i % 5)] = m.mle_expression(
                    &mut self.expressions,
                    &self.miscalls_exact,
                    &self.miscalls_mismatch,
                    &self.feature_models
                ) as f64;
            }

            // calculate mean change per feature of last 5 steps
            let mean_changes = last_changes.mean(Axis(1));
            debug!("mean changes={:?}", mean_changes);
            debug!("misidenfication probes={:?}", self.feature_models.iter().filter_map(|m| {
                if !m.feature_type.is_expressed() {
                    Some((m.counts.exact, self.miscalls_exact.total_to(m.feature_id), m.counts.mismatch, self.miscalls_mismatch.total_to(m.feature_id)))
                } else { None }
            }).collect_vec());
            //debug!("mean mean change={}", mean_changes.mean(Axis(0))[(0)]);

            assert_eq!(self.miscalls_exact.total_to(self.noise_id), 0);
            assert_eq!(self.miscalls_mismatch.total_to(self.noise_id), 0);

            debug!(
                "noise rate={}, exact miscalls={}, mismatch miscalls={}",
                self.expressions[self.noise_id],
                self.miscalls_exact.total_from(self.noise_id),
                self.miscalls_mismatch.total_from(self.noise_id)
            );

            eprintln!("x={:?}", self.expressions);
            //eprintln!("x={:?}", self.expressions.slice(s![..10]));

            // let mut likelihood = LogProb::ln_one();
            // for m in &self.feature_models {
            //     let x = self.expressions[m.feature_id];
            //     eprintln!("{} {:?} {} {:?} {:?}", x, self.miscalls_exact.miscalls.row(m.feature_id), m.feature_id, m.feature_type, m.event_probs);
            //     let l = m.likelihood(x, &self.miscalls_exact, &self.miscalls_mismatch);
            //     likelihood = likelihood + l;
            // }
            //
            // if let Some(debug_id) = self.debug_id {
            //     debug!(
            //         "DEBUG feature: x={}, miscalls_exact_from={}, miscalls_mismatch_from={}",
            //         self.expressions[debug_id],
            //         self.miscalls_exact.miscalls.row(debug_id),
            //         self.miscalls_mismatch.miscalls.row(debug_id)
            //     );
            // }
            //
            // //mem::swap(&mut self.expressions, &mut last_expressions);
            // if likelihood > map_likelihood {
            //     map_miscalls_exact = Some(self.miscalls_exact.clone());
            //     map_miscalls_mismatch = Some(self.miscalls_mismatch.clone());
            //     map_expressions = Some(self.expressions.clone());
            // }

            debug!("EM-iteration {} of {}", i, n_iterations);

            i += 1;
        }
        // update miscalls according to expressions of last iteration

        // self.miscalls_exact = map_miscalls_exact.unwrap();
        // self.miscalls_mismatch = map_miscalls_mismatch.unwrap();
        // self.expressions = map_expressions.unwrap();
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
