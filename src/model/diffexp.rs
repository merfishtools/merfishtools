use bio::stats::logprobs::LogProb;

use model;

pub type DiffexpMeasure = f64;
pub type CDF = model::dist::CDF<DiffexpMeasure>;


/// An estimate of differential expression
pub struct Estimate {
    pub differential_expression_pep: LogProb,
    pub differential_expression_bf: f64,
    pub expected_value: DiffexpMeasure,
    pub standard_deviation: DiffexpMeasure,
    pub credible_interval: (DiffexpMeasure, DiffexpMeasure)
}


impl CDF {
    /// Posterior error probability for differential expression.
    pub fn differential_expression_pep(&self, max_null_value: DiffexpMeasure) -> LogProb {
        self.get(&max_null_value).unwrap()
    }

    /// Bayes factor for differential expression.
    pub fn differential_expression_bf(&self, max_null_value: DiffexpMeasure) -> model::BayesFactor {
        let m0 = self.get(&max_null_value).unwrap();
        let m1 = 1.0 - m0;

        (m1 - m0).exp()
    }

    pub fn estimate(&self, max_fc: DiffexpMeasure) -> Estimate {
        let (ci_lower, ci_upper) = self.credible_interval();
        Estimate {
            differential_expression_pep: self.differential_expression_pep(max_fc),
            differential_expression_bf: self.differential_expression_bf(max_fc),
            expected_value: self.expected_value(),
            standard_deviation: self.standard_deviation(),
            credible_interval: (*ci_lower, *ci_upper)
        }
    }
}
