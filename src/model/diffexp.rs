use itertools::Itertools;
use bio::stats::logprobs::LogProb;
use bio::stats::logprobs;

use model;

pub type DiffexpMeasure = f64;
pub type PMF = model::pmf::PMF<DiffexpMeasure>;


/// An estimate of differential expression
pub struct Estimate {
    pub differential_expression_pep: LogProb,
    pub differential_expression_bf: f64,
    pub expected_value: DiffexpMeasure,
    pub standard_deviation: DiffexpMeasure,
    pub credible_interval: (DiffexpMeasure, DiffexpMeasure)
}


impl PMF {
    /// Posterior error probability for differential expression.
    pub fn differential_expression_pep(&self, max_null_value: DiffexpMeasure) -> LogProb {
        let probs = self.iter().filter(|e| e.value.abs() > max_null_value).map(|e| e.prob).collect_vec();
        logprobs::ln_1m_exp(logprobs::sum(&probs))
    }

    /// Bayes factor for differential expression.
    pub fn differential_expression_bf(&self, max_null_value: DiffexpMeasure) -> model::BayesFactor {
        let m0 = self.iter().filter(|e| e.value.abs() <= max_null_value).map(|e| e.prob).collect_vec();
        let m1 = self.iter().filter(|e| e.value.abs() > max_null_value).map(|e| e.prob).collect_vec();

        (logprobs::sum(&m1) - logprobs::sum(&m0)).exp()
    }

    pub fn estimate(&self, max_fc: DiffexpMeasure) -> Estimate {
        Estimate {
            differential_expression_pep: self.differential_expression_pep(max_fc),
            differential_expression_bf: self.differential_expression_bf(max_fc),
            expected_value: self.expected_value(),
            standard_deviation: self.standard_deviation(),
            credible_interval: self.credible_interval()
        }
    }
}
