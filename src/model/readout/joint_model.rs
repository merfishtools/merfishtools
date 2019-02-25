use std::collections::HashMap;

use itertools::Itertools;
use ndarray::prelude::*;
use ordered_float::NotNaN;
use rand;
use rand::distributions::IndependentSample;
use rand::{Rng, SeedableRng};

use bio::stats::{LogProb, Prob};

use crate::io::codebook::{Codebook, FeatureID};

use crate::model::readout::feature_model::AbstractFeatureModel;
use crate::model::readout::{Counts, Expressions, FeatureModel, Miscalls, NoiseModel, Xi};

pub struct JointModel {
    feature_models: HashMap<FeatureID, FeatureModel>,
    noise_model: NoiseModel,
    expressions: Expressions,
    miscalls_exact: Miscalls,
    miscalls_mismatch: Miscalls,
    margin: u32,
    em_run: bool,
    rng: rand::StdRng,
}

impl JointModel {
    pub fn new<'a, I: Iterator<Item = (&'a String, &'a Counts)>>(
        counts: I,
        p0: &[Prob],
        p1: &[Prob],
        codebook: &Codebook,
        window_width: u32,
        seed: usize,
    ) -> Self {
        let xi = Xi::new(p0, p1);
        let mut rng = rand::StdRng::from_seed(&[seed]);

        // Generate feature models and info about not expressed features.
        let mut feature_models = HashMap::new();
        let mut not_expressed_counts = Vec::new();
        let mut not_expressed_feature_ids = Vec::new();
        for (feature, counts) in counts {
            let feature = codebook.get_id(feature);
            let rec = codebook.record(feature);
            if rec.expressed() {
                feature_models.insert(
                    feature,
                    FeatureModel::new(feature, counts.clone(), codebook, &xi),
                );
            } else {
                not_expressed_counts.push(counts.clone());
                not_expressed_feature_ids.push(feature);
            }
        }

        let feature_count = codebook.len() + 1;

        // calculate start values (we take random numbers as start expression)
        let start_range = rand::distributions::Range::new(1, 10000);
        let mut expressions =
            Array1::from_iter((0..feature_count).map(|_| start_range.ind_sample(&mut rng)));
        for &feat_id in &not_expressed_feature_ids {
            expressions[feat_id] = 0;
        }

        // Take value larger than maximum feature id as noise id.
        let noise_id = codebook.len();
        let noise_model = NoiseModel::new(
            noise_id,
            not_expressed_feature_ids,
            not_expressed_counts,
            codebook,
            &xi,
        );

        let miscalls_exact = Miscalls::new(feature_count, &feature_models, &noise_model, true);
        let miscalls_mismatch = Miscalls::new(feature_count, &feature_models, &noise_model, false);

        JointModel {
            feature_models,
            noise_model,
            expressions,
            miscalls_exact,
            miscalls_mismatch,
            margin: window_width / 2,
            em_run: false,
            rng,
        }
    }

    /// Jointly estimate MAP using the EM algorithm.
    ///
    /// Expressions are our model parameters to estimate, miscalls are our latent variables,
    /// counts are our observed variables.
    ///
    /// # Arguments
    ///
    /// * ID of cell that is investigated (for debugging purposes).
    pub fn expectation_maximization(&mut self, cell: &str) {
        let mut feature_models: Vec<Box<&AbstractFeatureModel>> = self
            .feature_models
            .values()
            .map(|m| Box::new(m as &AbstractFeatureModel))
            .collect_vec();
        feature_models.push(Box::new(&self.noise_model));

        let n_iterations = 100;
        let change_window = 10;

        let mut last_changes = Array2::from_elem((self.expressions.len(), change_window), 0.0);

        // EM iterations
        for i in 1..n_iterations + 1 {
            debug!("EM-iteration {} of at most {}", i, n_iterations);

            // E-step: estimate miscalls
            debug!("E-step");

            // Shuffle models, such that miscalls can be drawn unbiased
            // this is necessary because there is a maximum miscall count for each feature.
            // If we would not shuffle, it would be easier for the first model to provide miscalls.
            self.rng.shuffle(&mut feature_models);
            for m in &feature_models {
                m.mle_miscalls(
                    &self.expressions,
                    &mut self.miscalls_exact,
                    &mut self.miscalls_mismatch,
                    &mut self.rng,
                );
            }

            // M-step: estimate expressions that maximize the probability
            debug!("M-step");
            for m in &feature_models {
                let feature_id = m.feature_id();
                last_changes[(feature_id, feature_id % change_window)] = m.mle_expression(
                    &mut self.expressions,
                    &self.miscalls_exact,
                    &self.miscalls_mismatch,
                ) as f64;
            }

            // calculate mean change per feature of last 5 steps
            let mean_changes = last_changes.mean_axis(Axis(1));
            debug!("TYPE=EM-iteration, CELL={}, x={:?}", cell, self.expressions);
            debug!("mean changes={:?}", mean_changes);
            let convergence = *mean_changes
                .iter()
                .map(|&c| NotNaN::new(c.abs()).unwrap())
                .max()
                .unwrap() <= 1.0;

            if i >= change_window && convergence {
                debug!("Convergence reached, stopping EM algorithm.");
                break;
            }
        }
        self.em_run = true;
        debug!(
            "TYPE=EM-result, CELL={}, noise-rate={}",
            cell,
            self.noise_rate()
        );
    }

    pub fn noise_rate(&self) -> f64 {
        assert!(self.em_run);
        self.expressions[self.expressions.len() - 1] as f64
            / self.expressions.slice(s![..-1]).scalar_sum() as f64
    }

    pub fn likelihood(&mut self, feature_id: FeatureID, x: u32) -> LogProb {
        assert!(self.em_run);
        self.feature_model(feature_id)
            .likelihood(x, &self.miscalls_exact, &self.miscalls_mismatch)
    }

    pub fn map_estimate(&self, feature_id: FeatureID) -> u32 {
        assert!(self.em_run);
        self.expressions[feature_id as usize]
    }

    fn feature_model(&self, feature_id: FeatureID) -> &FeatureModel {
        &self.feature_models[&feature_id]
    }

    /// Prior window for calculating expression PMF
    pub fn window(&self, feature_id: FeatureID) -> (u32, u32) {
        assert!(self.em_run);
        let x = self.map_estimate(feature_id);
        let xmin = self
            .feature_model(feature_id)
            .min_expression(&self.miscalls_exact, &self.miscalls_mismatch);
        let xmax = x + self.margin;
        assert!(
            xmax >= xmin,
            "xmax has to be greater than xmin: xmin={}, xmax={}",
            xmin,
            xmax
        );
        (xmin, xmax)
    }
}
