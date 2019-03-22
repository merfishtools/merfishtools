// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio::stats::Prob;
use clap::AppSettings::{ColoredHelp, DeriveDisplayOrder};
use clap::arg_enum;
use clap::_clap_count_exprs;
use failure::Error;
use itertools::Itertools;
use ordered_float::NotNaN;
use regex::Regex;
use structopt::StructOpt;

use merfishtools::cli;
use merfishtools::codebook;
use merfishtools::io::merfishdata;
use merfishtools::cli::Expression;

#[derive(StructOpt)]
#[structopt(
name = "merfishtools",
raw(global_settings = "&[ColoredHelp, DeriveDisplayOrder]")
)]
struct Opt {
    #[structopt(long, short, name = "verbose", parse(from_occurrences))]
    /// Level of verbosity. Can be specified multiple times.
    verbosity: usize,
    #[structopt(subcommand)]
    subcommand: Command,
}

#[derive(StructOpt)]
#[structopt(rename_all = "kebab-case")]
enum Command {
    #[structopt(
    name = "exp",
    about = "Estimate expressions for each feature (e.g. gene or transcript) in each cell.",
    after_help = "\
This command estimates expressions for each feature (e.g. gene or transcript) in each cell.
Results are provided as PMF (probability mass function) in columns:

cell
feature (e.g. gene, rna)
expression
posterior probability

Example usage:

merfishtools exp codebook.txt < data.txt > expression.txt"
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

        /// Path to write expected value and standard deviation estimates of expression to.
        ///
        //  Output is formatted into columns: cell, feature, expected value, standard deviation
        #[structopt(value_name = "TSV-FILE")]
        estimate: Option<String>,

        /// Path to write global statistics per cell to.
        ///
        //  Output is formatted into columns: cell, noise-rate
        #[structopt(long, value_name = "TSV-FILE")]
        stats: Option<String>,

        /// Seed for shuffling that occurs in EM algorithm.
        #[structopt(long, value_name = "INT")]
        seed: u64,

        /// Prior probability of 0->1 error
        #[structopt(
        long,
        default_value = "0.04",
        value_name = "FLOAT",
        multiple = true,
        )]
        p0: Vec<f64>,

        /// Prior probability of 1->0 error
        #[structopt(
        long,
        default_value = "0.10",
        value_name = "FLOAT",
        multiple = true,
        )]
        p1: Vec<f64>,

        /// Regular expression to select cells from cell column (see above).
        #[structopt(long, default_value = ".*", value_name = "REGEX")]
        cells: String,

        /// Width of the window to calculate PMF for.
        #[structopt(long, default_value = "100", value_name = "INT")]
        pmf_window_width: u32,

        #[structopt(long, default_value = "Bayes", value_name = "MODE", case_insensitive = true, raw(possible_values = "&Mode::variants()"))]
        mode: Mode,

        /// Number of threads to use.
        #[structopt(long, short, default_value = "1", value_name = "INT")]
        threads: usize,
    },

    #[structopt(
    name = "diffexp",
    about = "Test for differential expression between two groups of cells.",
    after_help = "\
This command calculates, for given expression PMFs (generated with merfishtools exp), differentially expressed features (e.g. genes or transcripts) between groups of cells given as separate input data.
Results are provided as columns:

feature (e.g. gene, rna)
posterior error probability (PEP) for differential expression
expected FDR when selecting all features down to the current
bayes factor (BF) for differential expression
expected log2 fold change of first vs second group
standard deviation of log2 fold change
lower and upper bound of 95% credible interval of log2 fold change

Example usage:

merfishtools diffexp data1.txt data2.txt > diffexp.txt"
    )]
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

    #[structopt(
    name = "multidiffexp",
    about = "Test for differential expression between multiple groups of cells.",
    after_help = "\
This command calculates, for given expression PMFs (obtained with merfishtools exp), differentially expressed features (e.g. genes or transcripts) between groups of cells given as separate input data.
Results are provided as columns:

feature (e.g. gene, rna)
posterior error probability (PEP) for differential expression
expected FDR when selecting all features down to the current
bayes factor (BF) for differential expression
expected coefficient of variation (CV)
standard deviation of CV
lower and upper bound of 95% credible interval of CV

Example usage:

merfishtools multidiffexp data1.txt data2.txt data3.txt > diffexp.txt"
    )]
    MultiDifferentialExpression {
        #[structopt(multiple = true)]
        /// Paths to expression PMFs for groups of cells.
        groups: Vec<String>,

        #[structopt(long, value_name = "FLOAT", default_value = "0.5")]
        /// Maximum coefficient of variation (CV) considered as no differential expression
        max_null_cv: f64,

        #[structopt(long, default_value = "1", value_name = "FLOAT")]
        /// Pseudocounts to add to means before CV calculation.
        pseudocounts: f64,

        #[structopt(long, value_name = "FILE")]
        /// Path to write CDFs of CVs to.
        cdf: Option<String>,

        #[structopt(long, short, default_value = "1", value_name = "INT")]
        /// Number of threads to use.
        threads: usize,
    },
    #[structopt(
    name = "est-error-rates",
    about = "Estimate 0-1 and 1-0 error rates.",
    after_help = "\
This command estimates 0-1 and 1-0 error rates from given MERFISH
readouts.

Example usage:

merfishtools est-error-rates readouts.tsv > error-rates.tsv

The produced output will have the three columns

pos
p0
p1

representing the position in the binary word, the 0-1 error rate and
the 1-0 error rate."
    )]
    EstimateErrors {
        #[structopt(value_name = "TSV-FILE")]
        codebook: String,

        #[structopt(value_name = "RAW-DATA")]
        ///Raw data containing molecule assignments to positions.
        ///
        /// If given as TSV file (ending on .tsv), the following columns are expected:
        /// cell, feature, readout
        /// Otherwise, the official MERFISH binary format is expected.
        raw_data: String,
    },
    #[structopt(
    name = "gen-mhd4",
    about = "Generate MERFISH MHD4 codebook with given parameters.",
    after_help = "\
This command generates a codebook with the given parameters.
Currently, the number of bits (N) is fixed to 16.

Example usage:

merfishtools gen-mhd4 -m 8 < transcript-names.txt > codebook.tsv

The output file codebook.tsv will contain the columns

feature (e.g. gene or transcript)
codeword"
    )]
    GenMhd4 {
        #[structopt(long, short = "m", value_name = "INT", required = true)]
        /// Number of 1-bits.
        onebits: u8,

        #[structopt(long, value_name = "PATTERN")]
        /// Regular expression pattern for features that should be marked
        /// as not expressed. This is useful to correctly model, e.g.,
        /// misidentification probes.
        not_expressed: Option<String>,
    },
    #[structopt(
    name = "gen-mhd2",
    about = "Generate MERFISH MHD2 codebook with given parameters.",
    after_help = "\
This command generates a codebook with the given parameters.

Example usage:

merfishtools gen-mhd2 -m 8 -N 16 < transcript-names.txt > codebook.tsv

The output file codebook.tsv will contain the columns

feature (e.g. gene or transcript)
codeword"
    )]
    GenMhd2 {
        #[structopt(long, short = "N", value_name = "INT")]
        /// Number of bits.
        bits: u8,

        #[structopt(long, short = "m", value_name = "INT")]
        /// Number of 1-bits.
        onebits: u8,

        #[structopt(long, value_name = "PATTERN")]
        /// Regular expression pattern for features that should be marked
        /// as not expressed. This is useful to correctly model, e.g.,
        /// misidentification probes.
        not_expressed: Option<String>,
    },
}

arg_enum! {
    #[derive(Debug)]
    enum Mode {
        Bayes,
        LA,
    }
}

#[allow(non_snake_case)]
fn main() -> Result<(), Error> {
    let opt = Opt::from_args(); // .version(env!("CARGO_PKG_VERSION"));

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
            mode,
            threads,
        } => {
            let convert_err_rates = |values: Vec<f64>| match values.len() {
                1 => vec![Prob::checked(values[0]).unwrap(); 32],
                _ => values
                    .into_iter()
                    .map(|p| Prob::checked(p).unwrap())
                    .collect_vec(),
            };
            match mode {
                Mode::Bayes => {
                    let mut expression = cli::ExpressionJBuilder::default()
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
                    expression.infer()
                }
                Mode::LA => {
                    let mut expression = merfishtools::model::la::expression::ExpressionTBuilder::default()
                        .p0(convert_err_rates(p0))
                        .p1(convert_err_rates(p1))
                        .codebook_path(codebook.to_owned())
                        .estimate_path(estimate.map(|v| v.to_owned()))
                        .stats_path(stats.map(|v| v.to_owned()))
                        .threads(threads)
                        .cells(Regex::new(&cells)?)
                        .max_hamming_distance(4 as usize) // TODO introduce mhd option
                        .bits(16)
                        .mode(merfishtools::model::la::expression::Mode::ErrorsThenExpression)
                        .seed(seed)
                        .build()
                        .unwrap();
                    if let merfishdata::Format::Binary = merfishdata::Format::from_path(&raw_data) {
                        expression.load_counts(&mut merfishdata::binary::Reader::from_file(&raw_data)?)?;
                    } else {
                        expression.load_counts(&mut merfishdata::tsv::Reader::from_file(&raw_data)?)?;
                    }

                    expression.infer()
                }
            }
        }

        Command::DifferentialExpression {
            group1,
            group2,
            max_null_log2fc,
            pseudocounts,
            cdf,
            threads,
        } => {
            let max_fc =
                NotNaN::new(max_null_log2fc).expect("NaN not allowed for --max-null-log2fc.");
            let pmf_path = cdf.as_ref().map(String::as_str);
            cli::differential_expression(&group1, &group2, pmf_path, max_fc, pseudocounts, threads)
        }

        Command::MultiDifferentialExpression {
            groups,
            max_null_cv,
            pseudocounts,
            cdf,
            threads,
        } => {
            let max_cv = NotNaN::new(max_null_cv).expect("NaN not allowed for --max-null-cv.");
            let group_paths = &groups.iter().map(String::as_ref).collect_vec();
            let cdf_path = cdf.as_ref().map(String::as_str);
            cli::multi_differential_expression(group_paths, cdf_path, max_cv, pseudocounts, threads)
        }

        Command::EstimateErrors { codebook, raw_data } => {
            cli::estimate_error_rates(&raw_data, &codebook)
        }

        Command::GenMhd2 {
            bits,
            onebits,
            not_expressed,
        } => {
            let words = codebook::generate_mhd2(bits, onebits);
            let not_expressed_pattern = not_expressed.as_ref().map(String::as_ref);
            cli::gen_codebook(&words, not_expressed_pattern)
        }

        Command::GenMhd4 {
            onebits,
            not_expressed,
        } => {
            let words = codebook::generate_mhd4(onebits);
            let not_expressed_pattern = not_expressed.as_ref().map(String::as_ref);
            cli::gen_codebook(&words, not_expressed_pattern)
        }
    }
}
