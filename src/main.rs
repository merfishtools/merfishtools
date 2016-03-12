#![allow(non_snake_case)]

extern crate merfishtools;
extern crate argparse;


use std::io::{stdout, stderr};
use std::str::FromStr;
use argparse::{ArgumentParser, Store, List, StoreTrue, StoreOption};

#[macro_use]
extern crate log;
extern crate fern;
extern crate bio;
extern crate csv;
extern crate itertools;
extern crate num;
extern crate rustc_serialize;
extern crate simple_parallel;
extern crate crossbeam;
extern crate regex;
extern crate rgsl;

pub mod model;
pub mod io;
pub mod cli;


#[allow(non_camel_case_types)]
#[derive(Debug)]
enum Command {
    exp,
    diffexp,
    //stats,
    None
}


impl FromStr for Command {
    type Err = ();

    fn from_str(src: &str) -> Result<Command, ()> {
        return match src {
            "exp"     => Ok(Command::exp),
            "diffexp" => Ok(Command::diffexp),
            //"stats"   => Ok(Command::stats),
            _         => Err(()),
        };
    }
}


fn main() {
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

    let mut debug = false;
    let mut subcommand = Command::None;
    let mut args = vec![];
    {
        let mut ap = ArgumentParser::new();
        ap.set_description("MERFISHtools is an analysis toolkit for MERFISH data.");
        ap.refer(&mut debug).add_option(&["--debug"], StoreTrue, "Print debugging information.");
        ap.refer(&mut subcommand).required()
            .add_argument("command", Store, "Command to run (exp or diffexp)");
        ap.refer(&mut args)
            .add_argument("arguments", List, "Arguments for command");
        ap.stop_on_first_argument(true);
        ap.parse_args_or_exit();
    }

    if let Err(e) = fern::init_global_logger(
        logger_config,
        if debug { log::LogLevelFilter::Debug } else { log::LogLevelFilter::Info }
    ) {
        panic!("Failed to initialize global logger: {}", e);
    }

    args.insert(0, format!("subcommand {:?}", subcommand));
    match subcommand {
        Command::exp     => exp(args),
        Command::diffexp => diffexp(args),
        //Command::stats   => stats(args),
        Command::None    => {
            error!("Unknown subcommand.");
            std::process::exit(1);
        }
    }
}


fn exp(args: Vec<String>) {
    let mut N = 16;
    let mut m = 4;
    let mut p0 = 0.04;
    let mut p1 = 0.1;
    let mut dist = 4;
    let mut threads = 1;
    let mut estimate_path = None;
    let mut codebook_path: Option<String> = None;
    let mut cells = ".*".to_owned();
    let mut window_width = 100;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"For given MERFISH data, calculate expressions for each feature (e.g. gene) in each cell.
Results are provided as PMF (probability mass function) in columns:

    cell
    feature (e.g. gene, rna)
    expression
    posterior probability

Example: 'merfishtools exp < data.txt > expression.txt'"#
        );

        ap.refer(&mut codebook_path).add_option(&["--codebook", "-c"], StoreOption, r#"
Path to codebook definition consisting of tab separated columns: feature, codeword, expressed.
The last column denotes if a codeword is assigned to e.g. a gene for which expression can be expected.
Unless you have misidentification probes (see Chen et al. Science 2015), you will have only ones in this column.
"#).required();
        ap.refer(&mut estimate_path).add_option(&["--estimate"], StoreOption, r#"
Path to write expected value and standard deviation estimates of expression to.
Output is formatted into columns: cell, feature, expected value, standard deviation
"#);
        ap.refer(&mut N).add_option(&["-N"], Store, "Number of bits in readout, i.e., number of hybridization rounds (default: 16).");
        ap.refer(&mut m).add_option(&["-m"], Store, "Number of 1-bits in readout (default: 4).");
        ap.refer(&mut p0).add_option(&["--p0"], Store, "Prior probability of 0->1 error (default: 0.04).");
        ap.refer(&mut p1).add_option(&["--p1"], Store, "Prior probability of 1->0 error (default: 0.1).");
        ap.refer(&mut dist).add_option(&["--hamming-dist", "--dist"], Store, "Hamming distance between encodings (default: 4).");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        ap.refer(&mut cells).add_option(&["--cells"], Store, "Regular expression for cells to select (default: all).");
        ap.refer(&mut window_width).add_option(&["--pmf-window-width"], Store, "Width of the window to calculate PMF for (default: 100).");
        parse_args_or_exit(&ap, args);
    }
    cli::expression(N, m, p0, p1, dist, &codebook_path.unwrap(), estimate_path, threads, &cells, window_width);
}


fn diffexp(args: Vec<String>) {
    let mut threads = 1;
    let mut group1_path = "".to_owned();
    let mut group2_path = "".to_owned();
    let mut pmf_path: Option<String> = None;
    let mut min_fc = 1.0f64.log2();

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"For given MERFISH data, calculate differentially expressed features (e.g. genes) between groups of cells given as separate input data.
Results are provided as columns:

    feature (e.g. gene, rna)
    posterior error probability (PEP) for differential expression
    expected log2 fold change
    standard deviation of log2 fold change

Example: "merfishtools diffexp data1.txt data2.txt > diffexp.txt""#
        );

        ap.refer(&mut pmf_path).add_option(&["--pmf"], StoreOption,
r#"Path to write PMF (probability mass function) of Log2 fold change to.
Output is formatted into columns: feature, foldchange, posterior probability"#);
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use (default: 1).");
        ap.refer(&mut min_fc)
          .add_option(&["--min-log2fc"], Store, "Minimum absolute log2 fold change considered as differential expression (default: log2(1.5)).");
        ap.refer(&mut group1_path).required()
          .add_argument("group1", Store, "Path to expression PMFs for group of cells.");
        ap.refer(&mut group2_path).required()
          .add_argument("group2", Store, "Path to expression PMFs for group of cells.");
        parse_args_or_exit(&ap, args);

    }
    cli::differential_expression(&group1_path, &group2_path, pmf_path, min_fc, threads);
}


/*fn stats(args: Vec<String>) {
    let mut N = 16;
    let mut m = 4;
    let mut p0 = 0.04;
    let mut p1 = 0.1;
    let mut dist = 4;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"For given MERFISH data, calculate per cell statistics.
Results are provided in columns:

cell
expected total expression

Example: 'merfishtools stats < stats.txt > expression.txt'"#
        );

        ap.refer(&mut N).add_option(&["-N"], Store, "Number of bits in readout, i.e., number of hybridization rounds (default: 16).");
        ap.refer(&mut m).add_option(&["-m"], Store, "Number of 1-bits in readout (default: 4).");
        ap.refer(&mut p0).add_option(&["--p0"], Store, "Prior probability of 0->1 error (default: 0.04).");
        ap.refer(&mut p1).add_option(&["--p1"], Store, "Prior probability of 1->0 error (default: 0.1).");
        ap.refer(&mut dist).add_option(&["--hamming-dist", "--dist"], Store, "Hamming distance between encodings (default: 4).");
        parse_args_or_exit(&ap, args);
    }
    cli::stats(N, m, p0, p1, dist);
}*/


fn parse_args_or_exit(ap: &ArgumentParser, args: Vec<String>) {
    match ap.parse(args, &mut stdout(), &mut stderr()) {
        Ok(()) =>  {}
        Err(x) => {
            std::process::exit(x);
        }
    }
}
