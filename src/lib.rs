// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[macro_use]
extern crate approx;
#[macro_use]
extern crate derive_builder;
#[macro_use]
extern crate failure_derive;
#[macro_use]
extern crate log;
#[macro_use(s)]
extern crate ndarray;
#[macro_use]
extern crate serde_derive;
#[macro_use]
extern crate serde_with;

pub mod cli;
pub mod codebook;
pub mod error_rates;
pub mod io;
pub mod model;
pub mod simulation;