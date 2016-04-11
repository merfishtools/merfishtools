use itertools::Itertools;

use model;


#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MeanVar {
    pub mean: f64,
    pub var: f64
}


impl MeanVar {

    pub fn standard_deviation(&self) -> f64 {
        self.var.sqrt()
    }

    pub fn coefficient_of_variation(&self) -> f64 {
        self.standard_deviation() / self.mean
    }
}


pub fn cdf<T: PartialOrd, F: Fn(f64, f64) -> T>(cdfs: &[model::dist::CDF<f64>], value: F) -> model::dist::CDF<T> {
    let n = cdfs.len() as f64;
    let mut curr = Vec::new();
    let mut prev = {
        let mut pmf = Vec::new();
        for (&value, prob) in cdfs[0].iter_pmf() {
            pmf.push(((value, 0.0), prob));
        }
        model::dist::CDF::from_pmf(pmf)
    };

    for (k, cdf) in cdfs.iter().enumerate().skip(1) {
        debug!("Iteration {}", k);
        let k = k as f64 + 1.0;

        curr = Vec::new();
        for (&(m, s), p) in prev.iter_pmf() {
            for (&value, prob) in cdf.iter_pmf() {
                let p = p + prob;
                let mk = m + (value - m) / k;
                let sk = s + (value - m) * (value - mk);
                curr.push(((mk, sk), p));
            }
        }
        debug!("PMF len={}", curr.len());
        prev = model::dist::CDF::from_pmf(curr).sample(1000);
    }
    let pmf = prev.iter_pmf().map(|(&(m, s), p)| {
        (value(m, s / (n - 1.0)), p)
    }).collect_vec();
    model::dist::CDF::from_pmf(pmf)
}
