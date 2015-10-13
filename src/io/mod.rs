pub mod merfishdata;
pub mod expression;
pub mod foldchange;


use bio::stats::logprobs::LogProb;


#[derive(RustcEncodable, RustcDecodable)]
pub struct PMF<T> {
    pub value: T,
    pub prob: LogProb
}
