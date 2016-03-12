pub use self::readout::Readout;

pub mod readout;
pub mod expression;
pub mod expressionset;
pub mod foldchange;
pub mod pmf;

pub const MIN_PROB: f64 = -13.815510557964274; // = 0.000001f64.ln();
//pub const MIN_PROB: f64 = -32.23619130191664;
