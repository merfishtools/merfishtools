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
extern crate statrs;
#[macro_use]
extern crate approx;
#[macro_use]
extern crate clap;
#[macro_use(s)]
extern crate ndarray;
extern crate bit_vec;
extern crate bit_set;
extern crate petgraph;
extern crate ordered_float;

use std::process;

use clap::App;
use itertools::Itertools;
use ordered_float::NotNaN;

use bio::stats::Prob;

pub mod model;
pub mod io;
pub mod cli;
pub mod codebook;


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
        let p0 = values_t!(matches, "p0", f64).unwrap_or_else(|e| e.exit());
        let p1 = values_t!(matches, "p1", f64).unwrap_or_else(|e| e.exit());
        let estimate_path = matches.value_of("estimate");
        let codebook_path = matches.value_of("codebook").unwrap();
        let cells = matches.value_of("cells").unwrap();
        let window_width = value_t!(matches, "pmf-window-width", u32).unwrap_or_else(|e| e.exit());
        let threads = value_t!(matches, "threads", usize).unwrap_or_else(|e| e.exit());

        let convert_err_rates = |values: Vec<f64>| {
            if values.len() == 1 {
                vec![Prob::checked(values[0]).unwrap(); 32]
            } else {
                values.into_iter().map(|p| {
                    Prob::checked(p).unwrap()
                }).collect_vec()
            }
        };

        cli::expression(
            convert_err_rates(p0),
            convert_err_rates(p1),
            &codebook_path,
            estimate_path,
            threads,
            &cells,
            window_width
        );
    } else if let Some(matches) = matches.subcommand_matches("diffexp") {
        let group1_path = matches.value_of("group1").unwrap();
        let group2_path = matches.value_of("group2").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_fc = NotNaN::new(value_t!(matches, "max-null-log2fc", f64).unwrap_or(1.0)).expect("NaN not allowed for --max-null-log2fc.");
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::differential_expression(&group1_path, &group2_path, cdf_path, max_fc, pseudocounts, threads);
    } else if let Some(matches) = matches.subcommand_matches("multidiffexp") {
        let group_paths = matches.values_of("groups").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_cv = NotNaN::new(value_t!(matches, "max-null-cv", f64).unwrap_or(0.5)).expect("NaN not allowed for --max-null-cv.");
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::multi_differential_expression(&group_paths.collect_vec(), cdf_path, max_cv, pseudocounts, threads);
    } else if let Some(matches) = matches.subcommand_matches("gen-mhd4") {
        let m = value_t!(matches, "onebits", u8).unwrap();
        let words = codebook::generate_mhd4(m);
        if let Err(e) = cli::gen_codebook(&words) {
            error!("{}", e);
            process::exit(1);
        }
    } else if let Some(matches) = matches.subcommand_matches("gen-mhd2") {
        let n = value_t!(matches, "bits", u8).unwrap();
        let m = value_t!(matches, "onebits", u8).unwrap();
        let words = codebook::generate_mhd2(n, m);
        if let Err(e) = cli::gen_codebook(&words) {
            error!("{}", e);
            process::exit(1);
        }
    }
}
