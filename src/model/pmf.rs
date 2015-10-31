use std::collections::HashMap;
use std::collections::hash_map;
use std::iter;
use std::slice;

use num::traits::{cast, NumCast};

use bio::stats::logprobs::LogProb;


pub struct PMF<T: Clone + Copy> {
    inner: Vec<(T, LogProb)>
}


impl<T: Clone + Copy> PMF<T> {

    pub fn new(inner: Vec<(T, LogProb)>) -> Self {
        PMF { inner: inner }
    }

    pub fn iter(&self) -> slice::Iter<(T, LogProb)> {
        self.inner.iter()
    }

    /// Return maximum a posteriori probability estimate (MAP).
    pub fn map(&self) -> T {
        let (mut max_x, mut max_prob) = self.iter().cloned().next().unwrap();
        for &(x, prob) in self.iter() {
            if prob >= max_prob {
                max_x = x;
                max_prob = prob;
            }
        }
        max_x
    }
}


impl<T: NumCast + Clone + Copy> PMF<T> {
    pub fn expected_value(&self) -> f64 {
        self.iter().map(|&(x, prob)| cast::<T, f64>(x).unwrap() * prob.exp()).fold(0.0f64, |s, e| s + e)
    }

    pub fn variance(&self) -> f64 {
        let e = self.expected_value();
        self.iter().map(|&(x, prob)| (cast::<T, f64>(x).unwrap() - e).powi(2) * prob.exp()).fold(0.0, |s, e| s + e)
    }
}
