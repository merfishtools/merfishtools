use std::f64;
use std::iter;
use std::slice;

use itertools::Itertools;

use bio::stats::logprobs::{self, LogProb};


#[derive(Debug, Clone)]
pub struct CDF<T: PartialOrd> {
    inner: Vec<(T, LogProb)>
}


impl<T: PartialOrd> CDF<T> {
    pub fn new(mut entries: Vec<(T, LogProb)>) -> Self {
        entries.sort_by(|&(ref a, _), &(ref b, _)| a.partial_cmp(b).unwrap());
        for i in 1..entries.len() {
            if entries[i].0 == entries[i - 1].0 {
                entries[i - 1].1 = logprobs::add(entries[i - 1].1, entries[i].1);
                entries.remove(i);
            }
            else {
                entries[i].1 = logprobs::add(entries[i - 1].1, entries[i].1);
            }
        }
        let mut cdf = CDF {
            inner: entries
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

    pub fn get(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.0.partial_cmp(value).unwrap()) {
                Ok(i) => self.inner[i].1,
                Err(i) => if i > 0 { self.inner[i - 1].1 } else { f64::NEG_INFINITY }
            })
        }
    }

    pub fn get_pmf(&self, value: &T) -> Option<LogProb> {
        if self.inner.is_empty() {
            None
        }
        else {
            Some(match self.inner.binary_search_by(|e| e.0.partial_cmp(value).unwrap()) {
                Ok(i) => if i > 0 { logprobs::sub(self.inner[i].1, self.inner[i - 1].1) } else { self.inner[0].1 },
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
    fn test_cdf() {
        let mut pmf = Vec::new();
        for i in 0..10 {
            pmf.push((i, 0.1f64.ln()));
        }

        let cdf = CDF::new(pmf.clone());
        for (value, prob) in pmf {
            assert_ulps_eq!(prob, cdf.get_pmf(&value).unwrap(), epsilon = 0.0000000000001);
        }
        assert_relative_eq!(cdf.total_prob(), 1.0f64.ln());
        assert_relative_eq!(cdf.get(&1).unwrap(), 0.2f64.ln());
    }
}
