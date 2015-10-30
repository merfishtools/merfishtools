use std::collections::HashMap;
use std::collections::hash_map;
use std::iter;

use num::traits::{cast, NumCast};

use bio::stats::logprobs::LogProb;


/*pub struct PMF<T: NumCast + Clone + Copy> {
    inner: Vec<(T, LogProb)>
}


impl<T: NumCast + Clone + Copy> PMF<T> {

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

    pub fn expected_value(&self) -> f64 {
        self.iter().map(|&(x, prob)| cast::<T, f64>(x).unwrap() * prob.exp()).fold(0.0f64, |s, e| s + e)
    }

    pub fn variance(&self) -> f64 {
        let e = self.expected_value();
        self.iter().map(|&(x, prob)| (cast::<T, f64>(x).unwrap() - e).powi(2) * prob.exp()).fold(0.0, |s, e| s + e)
    }
}*/


pub trait PMF<'a, T: Clone + Copy + Sized, D: Iterator<Item=T> + Sized> {
    fn domain(&'a self) -> D;

    fn cast(value: T) -> f64;

    fn posterior_prob(&'a self, value: T) -> LogProb;

    /// Return maximum a posteriori probability estimate (MAP).
    fn map(&'a self) -> T {
        let mut max_x = self.domain().next().unwrap();
        let mut max_prob = self.posterior_prob(max_x);
        for x in self.domain() {
            let prob = self.posterior_prob(x);
            if prob >= max_prob {
                max_x = x;
                max_prob = prob;
            }
        }
        max_x
    }

    fn expected_value(&'a self) -> f64 {
        self.domain().map(|x| Self::cast(x) * self.posterior_prob(x).exp()).fold(0.0f64, |s, e| s + e)
    }

    fn variance(&'a self) -> f64 {
        let e = self.expected_value();
        self.domain().map(|x| (Self::cast(x) - e).powi(2) * self.posterior_prob(x).exp()).fold(0.0, |s, e| s + e)
    }
}
