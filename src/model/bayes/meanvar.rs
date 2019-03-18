use ordered_float::NotNaN;

use bio::stats::probs::cdf;
use bio::stats::probs::cdf::CDF;
use itertools::Itertools;

#[derive(Debug, Clone, PartialOrd, PartialEq)]
pub struct MeanVar {
    pub mean: NotNaN<f64>,
    pub var: NotNaN<f64>,
}

impl MeanVar {
    pub fn standard_deviation(&self) -> f64 {
        self.var.sqrt()
    }

    pub fn coefficient_of_variation(&self) -> f64 {
        self.standard_deviation() / *self.mean
    }
}

pub fn cdf<T: Ord, F: Fn(NotNaN<f64>, NotNaN<f64>) -> T>(
    cdfs: &[CDF<NotNaN<f64>>],
    value: F,
) -> CDF<T> {
    assert!(
        cdfs.len() > 1,
        "meanvar::cdf() has to run on at least 2 cdfs"
    );
    let n = cdfs.len() as f64;
    let mut curr = Vec::new();
    let mut prev = {
        let mut pmf = Vec::new();
        for e in cdfs[0].iter_pmf() {
            pmf.push(cdf::Entry {
                value: (*e.value, NotNaN::new(0.0).unwrap()),
                prob: e.prob,
            });
        }
        CDF::from_pmf(pmf)
    };

    for (k, cdf) in cdfs.iter().enumerate().skip(1) {
        let k = k as f64 + 1.0;

        for prev_entry in prev.sample(1000).iter_pmf() {
            let &(m, s) = prev_entry.value;
            for e in cdf.iter_pmf() {
                let p = prev_entry.prob + e.prob;
                let mk = m + (*e.value - m) / k;
                let sk = s + (*e.value - m) * (*e.value - mk);
                curr.push(cdf::Entry {
                    value: (mk, sk),
                    prob: p,
                });
            }
        }
        prev = CDF::from_pmf(curr);
        curr = Vec::new();
    }
    let pmf = prev
        .iter_pmf()
        .map(|e| {
            let &(m, s) = e.value;
            cdf::Entry {
                value: value(m, s / (n - 1.0)),
                prob: e.prob,
            }
        })
        .collect_vec();
    CDF::from_pmf(pmf)
}

#[cfg(test)]
mod tests {
    use super::*;

    use bio::stats::LogProb;

    #[test]
    fn test_cdf() {
        let cdf1 = CDF::from_pmf(vec![
            cdf::Entry {
                value: NotNaN::new(6.0).unwrap(),
                prob: LogProb(-1.0),
            },
            cdf::Entry {
                value: NotNaN::new(70.0).unwrap(),
                prob: LogProb(-0.45867514538708193),
            },
        ]);
        let cdf2 = CDF::from_pmf(vec![
            cdf::Entry {
                value: NotNaN::new(3.0).unwrap(),
                prob: LogProb(-1.0),
            },
            cdf::Entry {
                value: NotNaN::new(70.0).unwrap(),
                prob: LogProb(-0.45867514538708193),
            },
        ]);
        let cdf3 = CDF::from_pmf(vec![
            cdf::Entry {
                value: NotNaN::new(3.0).unwrap(),
                prob: LogProb(-1.0),
            },
            cdf::Entry {
                value: NotNaN::new(70.0).unwrap(),
                prob: LogProb(-0.45867514538708193),
            },
        ]);
        println!("{}", cdf1.map().unwrap());
        println!("{}", cdf2.map().unwrap());

        let cdf = cdf(&[cdf1, cdf2, cdf3], |mean, _| mean);

        println!("{:?}", cdf);
        println!("{}", cdf.map().unwrap());
        assert_eq!(**cdf.map().unwrap(), 70.0);
    }
}
