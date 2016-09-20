// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

extern crate merfishtools;
#[macro_use]
extern crate log;
extern crate fern;
extern crate bio;
extern crate csv;
#[macro_use]
extern crate itertools;
extern crate num;
extern crate rustc_serialize;
extern crate cue;
extern crate regex;
extern crate rgsl;
extern crate ord_subset;
#[macro_use]
extern crate approx;
#[macro_use]
extern crate clap;

use clap::App;
use itertools::Itertools;

pub mod model;
pub mod io;
pub mod cli;


#[allow(non_snake_case)]
fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = App::from_yaml(yaml)
                      .version(env!("CARGO_PKG_VERSION"))
                      .get_matches();

    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, level: &log::LogLevel, _: &log::LogLocation| {
            match level {
                &log::LogLevel::Debug => format!("DEBUG: {}", msg),
                _ => msg.to_owned()
            }
        }),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Debug,
    };

    if let Err(e) = fern::init_global_logger(
        logger_config,
        if matches.is_present("verbose") { log::LogLevelFilter::Debug } else { log::LogLevelFilter::Info }
    ) {
        panic!("Failed to initialize logger: {}", e);
    }

    if let Some(matches) = matches.subcommand_matches("exp") {
        let N = value_t!(matches, "N", u8).unwrap_or(16);
        let m = value_t!(matches, "m", u8).unwrap_or(4);
        let p0 = value_t!(matches, "p0", f64).unwrap_or(0.04);
        let p1 = value_t!(matches, "p1", f64).unwrap_or(0.1);
        let dist = value_t!(matches, "hamming-dist", u8).unwrap_or(4);
        let estimate_path = matches.value_of("estimate");
        let codebook_path = matches.value_of("codebook").unwrap();
        let cells = matches.value_of("cells").unwrap_or(".*");
        let window_width = value_t!(matches, "pmf-window-width", u32).unwrap_or(100);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::expression(N, m, p0, p1, dist, &codebook_path, estimate_path, threads, &cells, window_width);
    } else if let Some(matches) = matches.subcommand_matches("diffexp") {
        let group1_path = matches.value_of("group1").unwrap();
        let group2_path = matches.value_of("group2").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_fc = value_t!(matches, "max-null-log2fc", f64).unwrap_or(1.0);
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::differential_expression(&group1_path, &group2_path, cdf_path, max_fc, pseudocounts, threads);
    } else if let Some(matches) = matches.subcommand_matches("multidiffexp") {
        let group_paths = matches.values_of("groups").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_cv = value_t!(matches, "max-null-cv", f64).unwrap_or(0.5);
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::multi_differential_expression(&group_paths.collect_vec(), cdf_path, max_cv, pseudocounts, threads);
    }
}
