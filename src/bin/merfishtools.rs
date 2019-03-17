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
    #[structopt(
    name = "exp",
    about = "Estimate expressions for each feature (e.g. gene or transcript) in each cell.",
    after_help =
    "This command estimates expressions for each feature (e.g. gene or transcript) in each cell.
Results are provided as PMF (probability mass function) in columns:

cell
feature (e.g. gene, rna)
expression
posterior probability

Example usage:

merfishtools exp codebook.txt < data.txt > expression.txt
                ")]
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

        #[structopt(long, value_name = "TSV-FILE")]
        /// Path to write global statistics per cell to.
        ///
        //  Output is formatted into columns: cell, noise-rate
        stats: Option<String>,

        #[structopt(long, value_name = "INT")]
        /// Seed for shuffling that occurs in EM algorithm.
        seed: usize,

        #[structopt(long, default_value = "0.04", value_name = "FLOAT", raw(use_delimiter = "true"))]
        /// Prior probability of 0->1 error
        p0: Vec<f64>,

        #[structopt(long, default_value = "0.10", value_name = "FLOAT", raw(use_delimiter = "true"))]
        /// Prior probability of 1->0 error
        p1: Vec<f64>,

        #[structopt(long, default_value = ".*", value_name = "REGEX")]
        /// Regular expression to select cells from cell column (see above).
        cells: String,

        #[structopt(long, default_value = "100", value_name = "INT")]
        /// Width of the window to calculate PMF for.
        pmf_window_width: u32,

        #[structopt(long, short, default_value = "1", value_name = "INT")]
        /// Number of threads to use.
        threads: usize,
    },

    #[structopt(name = "diffexp", about = "Test for differential expression between two groups of cells.")]
    DifferentialExpression {
        /// Path to expression PMFs for group of cells.
        group1: String,

        /// Path to expression PMFs for group of cells.
        group2: String,

        #[structopt(long, default_value = "1", value_name = "FLOAT")]
        /// Maximum absolute log2 fold change considered as no differential expression.
        max_null_log2fc: f64,

        #[structopt(long, default_value = "1", value_name = "FLOAT")]
        /// Pseudocounts to add to means before fold change calculation.
        pseudocounts: f64,

        #[structopt(long, value_name = "FILE")]
        /// Path to write CDFs of log2 fold changes to.
        cdf: Option<String>,

        #[structopt(long, short, default_value = "1", value_name = "INT")]
        /// Number of threads to use.
        threads: usize,
    },
    #[structopt(name = "multidiffexp")]
    MultiDifferentialExpression {},
    #[structopt(name = "est-error-rates")]
    EstimateErrors {},
}

#[allow(non_snake_case)]
fn main() -> Result<(), Error> {
    let opt = Opt::from_args();  // .version(env!("CARGO_PKG_VERSION"));

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
        match opt.verbosity {
            0 => log::LogLevelFilter::Info,
            _ => log::LogLevelFilter::Debug,
        },
    ) {
        panic!("Failed to initialize logger: {}", e);
    }
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
        } => {
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
                .codebook_path(codebook.to_owned())
                .estimate_path(estimate.map(|v| v.to_owned()))
                .stats_path(stats.map(|v| v.to_owned()))
                .threads(threads)
                .cells(Regex::new(&cells)?)
                .window_width(pmf_window_width)
                .seed(seed)
                .build()
                .unwrap();
            if let merfishdata::Format::Binary = merfishdata::Format::from_path(&raw_data) {
                expression.load_counts(&mut merfishdata::binary::Reader::from_file(&raw_data)?)?;
            } else {
                expression.load_counts(&mut merfishdata::tsv::Reader::from_file(&raw_data)?)?;
            }

            return expression.infer();
        }
        Command::DifferentialExpression {
            group1,
            group2,
            max_null_log2fc,
            pseudocounts,
            cdf,
            threads, } => {
            let max_fc = NotNaN::new(max_null_log2fc)
                .expect("NaN not allowed for --max-null-log2fc.");
            let pmf_path: Option<&str> = cdf.as_ref().map(String::as_str);
            return cli::differential_expression(
                &group1,
                &group2,
                pmf_path,
                max_fc,
                pseudocounts,
                threads,
            );
        }
        _ => unimplemented!()
    }
    process::exit(0);
    let matches = App::from_yaml(load_yaml!("../cli.yaml")).get_matches();
    if let Some(matches) = matches.subcommand_matches("diffexp") {
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
