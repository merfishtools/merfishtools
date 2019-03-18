pub mod cv;
pub mod diffexp;
pub mod expression;
pub mod expressionset;
pub mod foldchange;
pub mod meanvar;
pub mod pmf;
pub mod readout;

use bio::stats::LogProb;

pub const MIN_PROB: LogProb = LogProb(-13.815510557964274); // = 0.000001f64.ln();

pub type BayesFactor = f64;
