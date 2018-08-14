use std::slice;

use itertools::Itertools;

use bio::stats::LogProb;

#[derive(Clone, Debug)]
pub struct Entry<T: Clone> {
    pub value: T,
    pub prob: LogProb,
}

#[derive(Clone, Debug)]
pub struct PMF<T: Clone> {
    inner: Vec<Entry<T>>,
}

impl<T: Clone + Sized> PMF<T> {
    /// Create a new PMF from sorted vector.
    pub fn new(inner: Vec<Entry<T>>) -> Self {
        PMF { inner }
    }

    pub fn iter(&self) -> slice::Iter<Entry<T>> {
        self.inner.iter()
    }

    pub fn cdf(&self) -> Vec<LogProb> {
        LogProb::ln_cumsum_exp(self.inner.iter().map(|e| e.prob)).collect_vec()
    }

    /// Return maximum a posteriori probability estimate (MAP).
    pub fn map(&self) -> &T {
        let mut max = self.iter().next().unwrap();
        for e in self.iter() {
            if e.prob >= max.prob {
                max = e;
            }
        }
        &max.value
    }
}

impl<T: Clone + Sized + Copy> PMF<T> {
    /// Return the 95% credible interval.
    pub fn credible_interval(&self) -> (T, T) {
        let cdf = self.cdf();
        let lower = cdf
            .binary_search_by(|p| p.partial_cmp(&LogProb(0.025f64.ln())).unwrap())
            .unwrap_or_else(|i| i);
        let upper = cdf
            .binary_search_by(|p| p.partial_cmp(&LogProb(0.975f64.ln())).unwrap())
            .unwrap_or_else(|i| i);

        (self.inner[lower].value, self.inner[upper].value)
    }
}

impl PMF<f32> {
    pub fn scale(&mut self, scale: f32) {
        for e in &mut self.inner {
            e.value *= scale;
        }
    }
}
