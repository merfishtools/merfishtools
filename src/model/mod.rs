pub use self::readout::Readout;

pub mod readout;
pub mod expression;
pub mod expressionset;
pub mod diffexp;
pub mod foldchange;
pub mod cv;
pub mod pmf;

pub const MIN_PROB: f64 = -13.815510557964274; // = 0.000001f64.ln();

pub type BayesFactor = f64;
