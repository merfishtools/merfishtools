// Copyright 2016 Johannes Köster, 2019 Till Hartmann.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use clap::AppSettings::{ColoredHelp, DeriveDisplayOrder};
use failure::Error;
use rayon::ThreadPoolBuilder;
use regex::Regex;
use structopt::StructOpt;

use merfishtools::cli::Expression;
use merfishtools::io::merfishdata;

#[derive(StructOpt)]
#[structopt(
name = "barnacleboy",
raw(global_settings = "&[ColoredHelp, DeriveDisplayOrder]")
)]
struct Opt {
    /// Level of verbosity. Can be specified multiple times.
    #[structopt(long, short, name = "verbose", parse(from_occurrences))]
    verbosity: usize,

    #[structopt(subcommand)]
    subcommand: Command,
}

#[derive(StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Command {
    #[structopt(
    name = "exp",
    about = "Estimate expressions for each feature (e.g. gene or transcript) in each cell."
    )]
    Expression {
        /// Path to codebook definition consisting of tab separated columns: feature, codeword.
        ///
        /// Misidentification probes (see Chen et al. Science 2015) should not be contained in the codebook.
        #[structopt(value_name = "CODEBOOK-TSV")]
        codebook: String,

        /// Raw readout data containing molecule assignments to positions.
        ///
        /// If given as TSV file (ending on .tsv), the following columns are expected:
        /// cell, feature, hamming_dist, cell_position_x, cell_position_y, rna_position_x, rna_position_y.
        /// Otherwise, the official MERFISH binary format is expected.
        #[structopt(value_name = "READOUTS")]
        raw_data: String,

        /// Seed for shuffling that occurs in EM algorithm.
        #[structopt(long, value_name = "INT")]
        seed: u64,

        /// Number of bits used for barcodes.
        #[structopt(long, value_name = "INT", default_value = "0")]
        num_bits: usize,

        /// Prior probability of 1 ← 0 error
        #[structopt(long, default_value = "0.04", value_name = "FLOAT", multiple = true)]
        p0: Vec<f64>,

        /// Prior probability of 0 ← 1 error
        #[structopt(long, default_value = "0.09", value_name = "FLOAT", multiple = true)]
        p1: Vec<f64>,

        /// Regular expression to select cells from cell column (see above).
        #[structopt(long, default_value = ".*", value_name = "REGEX")]
        cells: String,

        /// Number of threads to use. If 0, uses `RAYON_NUM_THREADS' environment variable.
        #[structopt(long, short, default_value = "0", value_name = "INT")]
        threads: usize,

        /// Maximum hamming distance to consider.
        #[structopt(long, short = "d", value_name = "INT", default_value = "5")]
        max_hamming_distance: usize,

        // Relaxation parameter for successive-over-relaxation step. 0 < omega < 2.
        #[structopt(long, short = "w", value_name = "FLOAT", default_value = "1.25")]
        omega: f32,

        /// Mode of operation: Estimate errors or expression first?
        #[structopt(
        long,
        short = "m",
        default_value = "ErrorsThenExpression",
        raw(possible_values = "&merfishtools::model::la::expression::Mode::variants()"),
        case_insensitive = true
        )]
        mode: merfishtools::model::la::expression::Mode,

        /// Path to write expected value of expression to.
        ///
        /// Output is formatted into columns: cell, feature, expected value
        #[structopt(long, short = "o", value_name = "TSV-FILE")]
        estimate: Option<String>,

        /// Path to write estimated transition probabilities to.
        ///
        /// TODO: Format prone to change.
        #[structopt(long, short = "e", value_name = "TSV-FILE")]
        errors: Option<String>,
    },
}

#[allow(non_snake_case)]
fn main() -> Result<(), Error> {
    let opt = Opt::from_args(); // .version(env!("CARGO_PKG_VERSION"));

    fern::Dispatch::new()
        .format(|out, msg, record| match record.level() {
            log::Level::Debug => out.finish(format_args!("DEBUG: {}", msg)),
            _ => out.finish(msg.to_owned()),
        })
        .level(match opt.verbosity {
            0 => log::LevelFilter::Info,
            _ => log::LevelFilter::Debug,
        })
        .chain(std::io::stderr())
        .apply()?;

    match opt.subcommand {
        // it is not yet possible to use `exp @ Command::Expression { .. } => { do stuff with exp }`
        // see https://github.com/rust-lang/rfcs/pull/2593
        // and https://github.com/varkor/rfcs/blob/enum-variant-types/text/0000-enum-variant-types.md
        Command::Expression {
            codebook,
            raw_data,
            seed,
            num_bits,
            p0,
            p1,
            cells,
            mode,
            threads,
            max_hamming_distance,
            omega,
            estimate,
            errors,
        } => {
            let convert_err_rates = |values: Vec<f64>| match values.len() {
                1 => vec![values[0]; num_bits],
                _ => values,
            };
            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()?;
            let mut expression = merfishtools::model::la::expression::ExpressionT::new(
                convert_err_rates(p0),
                convert_err_rates(p1),
                codebook.to_owned(),
                estimate.map(|v| v.to_owned()),
                errors.map(|v| v.to_owned()),
                Regex::new(&cells)?,
                max_hamming_distance,
                omega,
                mode,
                seed,
                num_bits,
            );
            match merfishdata::Format::from_path(&raw_data) {
                merfishdata::Format::Binary => expression.load_counts(
                    &mut merfishdata::binary::Reader::from_file(&raw_data)?,
                    merfishdata::Format::Binary,
                )?,
                merfishdata::Format::TSV => expression.load_counts(
                    &mut merfishdata::tsv::Reader::from_file(&raw_data)?,
                    merfishdata::Format::TSV,
                )?,
                merfishdata::Format::Simulation => expression.load_counts(
                    &mut merfishdata::sim::Reader::from_file(&raw_data)?,
                    merfishdata::Format::Simulation,
                )?,
            }
            expression.infer()
        }
    }
}
