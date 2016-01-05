use std::slice;

use num::traits::{cast, NumCast};

use bio::stats::logprobs::LogProb;
use bio::stats::logprobs;


#[derive(Clone, Debug)]
pub struct Entry<T: Clone + Copy> {
    pub value: T,
    pub prob: LogProb
}


#[derive(Clone, Debug)]
pub struct PMF<T: Clone + Copy> {
    inner: Vec<Entry<T>>
}


impl<T: Clone + Copy + Sized> PMF<T> {

    /// Create a new PMF from sorted vector.
    pub fn new(inner: Vec<Entry<T>>) -> Self {
        PMF { inner: inner }
    }

    pub fn iter(&self) -> slice::Iter<Entry<T>> {
        self.inner.iter()
    }

    pub fn cdf(&self) -> Vec<LogProb> {
        logprobs::cumsum(self.inner.iter().map(|e| e.prob))
    }

    /// Return maximum a posteriori probability estimate (MAP).
    pub fn map(&self) -> T {
        let mut max = self.iter().next().unwrap();
        for e in self.iter() {
            if e.prob >= max.prob {
                max = e;
            }
        }
        max.value
    }

    /// Return the 5%-95% credible interval.
    pub fn credible_interval(&self) -> (T, T) {
        let cdf = self.cdf();
        let lower = cdf.binary_search_by(|p| p.partial_cmp(&0.025).unwrap()).unwrap_or_else(|i| i);
        let upper = cdf.binary_search_by(|p| p.partial_cmp(&0.975).unwrap()).unwrap_or_else(|i| i - 1);

        (self.inner[lower].value, self.inner[upper].value)
    }
}


impl<T: NumCast + Clone + Copy> PMF<T> {
    pub fn expected_value(&self) -> f64 {
        self.iter().map(|e| {
            cast::<T, f64>(e.value).unwrap() * e.prob.exp()
        }).fold(0.0f64, |s, e| s + e)
    }

    pub fn variance(&self) -> f64 {
        let ev = self.expected_value();
        self.iter().map(|e| {
                (cast::<T, f64>(e.value).unwrap() - ev).powi(2) * e.prob.exp()
        }).fold(0.0, |s, e| s + e)
    }

    pub fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }
}


impl PMF<f32> {
    pub fn scale(&mut self, scale: f32) {
        for e in self.inner.iter_mut() {
            e.value *= scale;
        }
    }
}
