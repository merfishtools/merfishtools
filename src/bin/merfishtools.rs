// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#[macro_use]
extern crate clap;

use bio::stats::Prob;
use clap::App;
use failure::Error;
use itertools::Itertools;
use ordered_float::NotNaN;
use regex::Regex;

use merfishtools::cli;
use merfishtools::codebook;
use merfishtools::io::merfishdata;
use structopt::StructOpt;
use clap::AppSettings::{ColoredHelp, DeriveDisplayOrder};
use std::process;

#[derive(StructOpt)]
#[structopt(name = "merfishtools",
raw(global_settings = "&[ColoredHelp, DeriveDisplayOrder]"))]
struct Opt {
    #[structopt(long, short, name = "verbose", parse(from_occurrences))]
    /// Level of verbosity. Can be specified multiple times.
    verbosity: usize,
    #[structopt(subcommand)]
    subcommand: Command,
}

#[derive(StructOpt)]
enum Command {
    #[structopt(name = "exp")]
    Expression {
        #[structopt(value_name = "CODEBOOK-TSV")]
        /// Path to codebook definition consisting of tab separated columns: feature, codeword.
        ///
        /// Misidentification probes (see Chen et al. Science 2015) should not be contained in the codebook.
        codebook: String,

        #[structopt(value_name = "READOUTS")]
        /// Raw readout data containing molecule assignments to positions.
        ///
        /// If given as TSV file (ending on .tsv), the following columns are expected:
        /// cell, feature, hamming_dist, cell_position_x, cell_position_y, rna_position_x, rna_position_y.
        /// Otherwise, the official MERFISH binary format is expected.
        raw_data: String,

        #[structopt(value_name = "TSV-FILE")]
        /// Path to write expected value and standard deviation estimates of expression to.
        ///
        //  Output is formatted into columns: cell, feature, expected value, standard deviation
        estimate: Option<String>,

        #[structopt(long)]
        /// Path to write global statistics per cell to.
        ///
        //  Output is formatted into columns: cell, noise-rate
        stats: Option<String>,

        #[structopt(long)]
        /// Seed for shuffling that occurs in EM algorithm.
        seed: Option<usize>,

        #[structopt(long, default_value = "0.04", raw(use_delimiter = "true"))]
        /// Prior probability of 0->1 error
        p0: Vec<f64>,

        #[structopt(long, default_value = "0.10", raw(use_delimiter = "true"))]
        /// Prior probability of 1->0 error
        p1: Vec<f64>,

        #[structopt(long, default_value = ".*", value_name = "REGEX")]
        /// Regular expression to select cells from cell column (see above).
        cells: String,

        #[structopt(long, default_value = "100")]
        /// Width of the window to calculate PMF for.
        pmf_window_width: u32,

        #[structopt(long, short, default_value = "1")]
        /// Number of threads to use.
        threads: usize,
    },
    #[structopt(name = "diffexp")]
    DifferentialExpression {},
    #[structopt(name = "multidiffexp")]
    MultiDifferentialExpression {},
    #[structopt(name = "est-error-rates")]
    EstimateErrors {},
}

#[allow(non_snake_case)]
fn main() -> Result<(), Error> {
    let opt = Opt::from_args();
    match opt.subcommand {
        // it is not yet possible to use `exp @ Command::Expression { .. } => { do stuff with exp }`
        // see https://github.com/rust-lang/rfcs/pull/2593
        // and https://github.com/varkor/rfcs/blob/enum-variant-types/text/0000-enum-variant-types.md
        Command::Expression {
            codebook,
            raw_data,
            estimate,
            stats,
            seed,
            p0,
            p1,
            cells,
            pmf_window_width,
            threads
        } => unimplemented!(),
        _ => unimplemented!()
    }
    process::exit(0);
    let yaml = load_yaml!("../cli.yaml");
    let matches = App::from_yaml(yaml)
        .version(env!("CARGO_PKG_VERSION"))
        .get_matches();

    let logger_config = fern::DispatchConfig {
        format: Box::new(
            |msg: &str, level: &log::LogLevel, _: &log::LogLocation| match *level {
                log::LogLevel::Debug => format!("DEBUG: {}", msg),
                _ => msg.to_owned(),
            },
        ),
        output: vec![fern::OutputConfig::stderr()],
        level: log::LogLevelFilter::Debug,
    };

    if let Err(e) = fern::init_global_logger(
        logger_config,
        if matches.is_present("verbose") {
            log::LogLevelFilter::Debug
        } else {
            log::LogLevelFilter::Info
        },
    ) {
        panic!("Failed to initialize logger: {}", e);
    }

    if let Some(matches) = matches.subcommand_matches("exp") {
        let p0 = values_t!(matches, "p0", f64).unwrap_or_else(|e| e.exit());
        let p1 = values_t!(matches, "p1", f64).unwrap_or_else(|e| e.exit());
        let estimate_path = matches.value_of("estimate");
        let stats_path = matches.value_of("stats");
        let codebook_path = matches.value_of("codebook").unwrap();
        let cells = matches.value_of("cells").unwrap();
        let window_width = value_t!(matches, "pmf-window-width", u32).unwrap_or_else(|e| e.exit());
        let threads = value_t!(matches, "threads", usize).unwrap_or_else(|e| e.exit());
        let seed = value_t!(matches, "seed", usize).unwrap_or_else(|e| e.exit());
        let raw_data = matches.value_of("raw_data").unwrap();

        let convert_err_rates = |values: Vec<f64>| {
            if values.len() == 1 {
                vec![Prob::checked(values[0]).unwrap(); 32]
            } else {
                values
                    .into_iter()
                    .map(|p| Prob::checked(p).unwrap())
                    .collect_vec()
            }
        };

        let mut expression = cli::ExpressionBuilder::default()
            .p0(convert_err_rates(p0))
            .p1(convert_err_rates(p1))
            .codebook_path(codebook_path.to_owned())
            .estimate_path(estimate_path.map(|v| v.to_owned()))
            .stats_path(stats_path.map(|v| v.to_owned()))
            .threads(threads)
            .cells(Regex::new(cells)?)
            .window_width(window_width)
            .seed(seed)
            .build()
            .unwrap();
        if let merfishdata::Format::Binary = merfishdata::Format::from_path(raw_data) {
            expression.load_counts(&mut merfishdata::binary::Reader::from_file(raw_data)?)?;
        } else {
            expression.load_counts(&mut merfishdata::tsv::Reader::from_file(raw_data)?)?;
        }

        expression.infer()
    } else if let Some(matches) = matches.subcommand_matches("diffexp") {
        let group1_path = matches.value_of("group1").unwrap();
        let group2_path = matches.value_of("group2").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_fc = NotNaN::new(value_t!(matches, "max-null-log2fc", f64).unwrap_or(1.0))
            .expect("NaN not allowed for --max-null-log2fc.");
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::differential_expression(
            &group1_path,
            &group2_path,
            cdf_path,
            max_fc,
            pseudocounts,
            threads,
        )
    } else if let Some(matches) = matches.subcommand_matches("multidiffexp") {
        let group_paths = matches.values_of("groups").unwrap();
        let cdf_path = matches.value_of("cdf");
        let max_cv = NotNaN::new(value_t!(matches, "max-null-cv", f64).unwrap_or(0.5))
            .expect("NaN not allowed for --max-null-cv.");
        let pseudocounts = value_t!(matches, "pseudocounts", f64).unwrap_or(1.0);
        let threads = value_t!(matches, "threads", usize).unwrap_or(1);

        cli::multi_differential_expression(
            &group_paths.collect_vec(),
            cdf_path,
            max_cv,
            pseudocounts,
            threads,
        )
    } else if let Some(matches) = matches.subcommand_matches("gen-mhd4") {
        let m = value_t!(matches, "onebits", u8).unwrap();
        let not_expressed_pattern = matches.value_of("not-expressed");
        let words = codebook::generate_mhd4(m);
        cli::gen_codebook(&words, not_expressed_pattern)
    } else if let Some(matches) = matches.subcommand_matches("gen-mhd2") {
        let n = value_t!(matches, "bits", u8).unwrap();
        let m = value_t!(matches, "onebits", u8).unwrap();
        let not_expressed_pattern = matches.value_of("not-expressed");
        let words = codebook::generate_mhd2(n, m);
        cli::gen_codebook(&words, not_expressed_pattern)
    } else if let Some(matches) = matches.subcommand_matches("est-error-rates") {
        let codebook = matches.value_of("codebook").unwrap();
        let raw_data = matches.value_of("raw_data").unwrap();

        cli::estimate_error_rates(raw_data, codebook)
    } else {
        panic!("bug: unexpected subcommand");
    }
}
