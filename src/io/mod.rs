// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

pub mod cdf;
pub mod codebook;
pub mod simple_codebook;
pub mod estimation;
pub mod merfishdata;

#[derive(RustcEncodable, RustcDecodable, PartialEq, Eq, Clone, Copy, Hash)]
pub struct Cell {
    pub experiment: u32,
    pub cell: u32,
}
