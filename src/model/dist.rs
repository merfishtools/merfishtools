use std::collections::{HashMap, hash_map};
use std::f64;
use std::hash::Hash;
use std::iter;
use std::slice;

use itertools::Itertools;

use bio::stats::logprobs::{self, LogProb};


#[derive(Debug, Clone)]
pub struct PMF<T: Eq + Hash> {
    inner: HashMap<T, LogProb>
}


impl<T: Eq + Hash> PMF<T> {
    pub fn new() -> Self {
        PMF { inner: HashMap::new() }
    }

    pub fn add(&mut self, value: T, prob: LogProb) {
        let p = self.inner.entry(value).or_insert(f64::NEG_INFINITY);
        *p = logprobs::add(*p, prob);
    }

    pub fn get(&self, value: &T) -> Option<LogProb> {
        self.inner.get(value).cloned()
    }

    pub fn clear(&mut self) {
        self.inner.clear()
    }

    pub fn iter(&self) -> hash_map::Iter<T, LogProb> {
        self.inner.iter()
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}


impl<T: PartialOrd + Eq + Hash> PMF<T> {
    pub fn to_cdf(self) -> CDF<T> {
        CDF::new(self)
    }
}


impl<T: Eq + Hash> IntoIterator for PMF<T> {
    type Item = (T, LogProb);
    type IntoIter = hash_map::IntoIter<T, LogProb>;

    fn into_iter(self) -> Self::IntoIter {
        self.inner.into_iter()
    }
}


#[derive(Debug, Clone)]
pub struct CDF<T: PartialOrd> {
    inner: Vec<(T, LogProb)>
}


impl<T: PartialOrd> CDF<T> {
    pub fn new<I: IntoIterator<Item = (T, LogProb)>>(entries: I) -> Self {
        let mut inner = entries.into_iter().sorted_by(|&(ref a, _), &(ref b, _)| a.partial_cmp(b).unwrap());
        {
            let mut prob_sum = f64::NEG_INFINITY;
            for e in inner.iter_mut() {
                e.1 = logprobs::add(e.1, prob_sum);
                prob_sum = e.1;
            }
        }
        let mut cdf = CDF {
            inner: inner
        };

        if relative_eq!(cdf.total_prob(), 0.0) && cdf.total_prob() > 0.0 {
            cdf.inner.last_mut().unwrap().1 = 0.0;
        }

        cdf
    }

    pub fn sample(mut self, n: usize) -> Self {
        if self.inner.len() <= n {
            println!("not sampling");
            self
        }
        else {
            let s = self.inner.len() / (n - 1);
            let last = self.inner.pop().unwrap();
            let mut inner = self.inner.into_iter().step(s).collect_vec();
            inner.push(last);
            CDF {
                inner: inner
            }
        }
    }

    pub fn get(&self, value: T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.0.partial_cmp(&value).unwrap()) {
                Ok(i) => self.inner[i].1,
                Err(i) => if i > 0 { self.inner[i - 1].1 } else { f64::NEG_INFINITY }
            })
        }
    }

    pub fn total_prob(&self) -> LogProb {
        self.inner.last().unwrap().1
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

impl<T: Clone + PartialOrd> CDF<T> {
    pub fn iter_pmf(&self) -> CDFPMFIter<T> {
        fn cdf_to_pmf<G: Clone>(last_prob: &mut LogProb, e: &(G, LogProb)) -> Option<(G, LogProb)> {
            let &(ref value, cdf_prob) = e;
            let prob = logprobs::sub(cdf_prob, *last_prob);
            *last_prob = cdf_prob;
            Some((value.clone(), prob))
        }
        self.inner.iter().scan(f64::NEG_INFINITY, cdf_to_pmf)
    }
}

pub type CDFPMFIter<'a, T> = iter::Scan<slice::Iter<'a, (T, LogProb)>, LogProb, fn(&mut LogProb, &(T, LogProb)) -> Option<(T, LogProb)>>;


#[cfg(test)]
mod test {
    use super::*;

    use itertools::Itertools;

    #[test]
    fn test_pmf() {
        let mut pmf = PMF::new();
        pmf.add(1, 0.5f64.ln());
        pmf.add(2, 0.1f64.ln());
        pmf.add(1, 0.25f64.ln());

        assert_relative_eq!(pmf.get(&1).unwrap(), 0.75f64.ln());
        assert_relative_eq!(pmf.get(&2).unwrap(), 0.1f64.ln());
    }

    #[test]
    fn test_cdf() {
        let mut pmf = PMF::new();
        for i in 0..10 {
            pmf.add(i, 0.1f64.ln());
        }

        let cdf = pmf.clone().to_cdf();
        for (value, prob) in cdf.iter_pmf() {
            assert_ulps_eq!(prob, pmf.get(&value).unwrap(), epsilon = 0.0000000000001);
        }
        assert_relative_eq!(cdf.total_prob(), 1.0f64.ln());
        assert_relative_eq!(cdf.get(1).unwrap(), 0.2f64.ln());
    }
}
