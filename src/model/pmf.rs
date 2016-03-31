use std::slice;

use itertools::Itertools;
use num::traits::{cast, NumCast, Zero, One};
use num::rational::Ratio;

use bio::stats::logprobs::LogProb;
use bio::stats::logprobs;

use model;


#[derive(Clone, Debug)]
pub struct Entry<T: Clone> {
    pub value: T,
    pub prob: LogProb
}


#[derive(Clone, Debug)]
pub struct PMF<T: Clone> {
    inner: Vec<Entry<T>>
}


impl<T: Clone + Sized> PMF<T> {

    /// Create a new PMF from sorted vector.
    pub fn new(inner: Vec<Entry<T>>) -> Self {
        PMF { inner: inner }
    }

    pub fn iter(&self) -> slice::Iter<Entry<T>> {
        self.inner.iter()
    }

    pub fn cdf(&self) -> Vec<LogProb> {
        logprobs::cumsum(self.inner.iter().map(|e| e.prob)).collect_vec()
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
    /// Return the 5%-95% credible interval.
    pub fn credible_interval(&self) -> (T, T) {
        let cdf = self.cdf();
        let lower = cdf.binary_search_by(|p| p.partial_cmp(&0.025f64.ln()).unwrap()).unwrap_or_else(|i| i);
        let upper = cdf.binary_search_by(|p| p.partial_cmp(&0.975f64.ln()).unwrap()).unwrap_or_else(|i| i - 1);

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


#[derive(Debug)]
pub struct MeanVar {
    pub mean: f64,
    pub var: f64,
    pub prob: f64
}


impl MeanVar {
    pub fn new(pmfs: &[PMF<Ratio<u32>>]) -> Vec<MeanVar> {
        let n = Ratio::from_integer(pmfs.len() as u32);
        let mut curr = model::dist::PMF::new();
        let mut prev = {
            let mut pmf = model::dist::PMF::new();
            for e in pmfs[0].iter() {
                pmf.add((e.value, Ratio::zero()), e.prob);
            }
            pmf.to_cdf()
        };



        for (k, pmf) in pmfs.iter().enumerate().skip(1) {
            debug!("Iteration {}", k);
            let k = Ratio::from_integer(k as u32 + 1);

            curr = model::dist::PMF::new();
            for ((m, s), p) in prev.iter_pmf() {
                for x in pmf.iter() {
                    let p = p + x.prob;
                    let mk = m + (x.value - m) / k;
                    let sk = s + (x.value - m) * (x.value - mk);
                    curr.add((mk, sk), p);
                }
            }
            debug!("PMF len={}", curr.len());
            prev = curr.to_cdf().sample(1000);
        }

        let to_f64 = |ratio: Ratio<u32>| *ratio.numer() as f64 / *ratio.denom() as f64;

        prev.iter_pmf().filter_map(|((m, s), p)| {
            if p >= model::MIN_PROB {
                Some(MeanVar { mean: to_f64(m), var: to_f64(s / (n - Ratio::one())), prob: p })
            }
            else {
                None
            }
        }).collect_vec()
    }

    pub fn standard_deviation(&self) -> f64 {
        self.var.sqrt()
    }

    pub fn coefficient_of_variation(&self) -> f64 {
        self.standard_deviation() / self.mean
    }
}
