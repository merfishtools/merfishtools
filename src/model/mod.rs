pub use self::readout::Readout;

pub mod readout;
pub mod expression;
pub mod expressionset;
pub mod foldchange;
pub mod pmf;

const MIN_PROB: f64 = -13.815510557964274; // = 0.000001f64.ln();
