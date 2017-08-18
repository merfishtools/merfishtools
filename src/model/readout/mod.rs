use std::cmp;
use std::collections::HashMap;

use ndarray::prelude::*;

use io::codebook::FeatureID;

pub mod joint_model;
pub mod feature_model;
pub mod xi;

pub use model::readout::joint_model::JointModel;
pub use model::readout::feature_model::{FeatureModel, NoiseModel};
pub use model::readout::xi::Xi;


pub type Expressions = Array1<u32>;



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
    pub fn new(
        feature_count: usize,
        feature_models: &HashMap<FeatureID, FeatureModel>,
        noise_model: &NoiseModel,
        exact: bool
    ) -> Self {
        let mut max_total_to = vec![0; feature_count];
        for m in feature_models.values() {
            max_total_to[m.feature_id] = if exact { m.counts.exact } else { m.counts.mismatch };
        }
        for (&feature_id, counts) in noise_model.not_expressed_feature_ids.iter().zip(
            &noise_model.not_expressed_counts
        ) {
            max_total_to[feature_id] = if exact { counts.exact } else { counts.mismatch };
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

    /// Set a given miscall value, such that the total maximum per feature capacity is not exceeded.
    /// The given capacity is updated.
    /// Returns the change.
    pub fn set(&mut self, i: FeatureID, j: FeatureID, value: u32, capacity: &mut u32) -> u32 {
        let change = {
            let entry = self.miscalls.get_mut((i as usize, j as usize)).unwrap();
            let value = cmp::min(
                cmp::min(value, *capacity),
                self.max_total_to[j] - (self.total_to[j] - *entry)
            );
            let change = value as i32 - *entry as i32;
            *capacity -= value;
            *entry = value;
            self.total_to[j] = (self.total_to[j] as i32 + change) as u32;
            self.total = (self.total as i32 + change) as u32;

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

    pub fn max_total_to(&self, feature: FeatureID) -> u32 {
        self.max_total_to[feature]
    }

    pub fn total(&self) -> u32 {
        self.total
    }
}
