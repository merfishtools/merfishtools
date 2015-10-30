#![allow(non_snake_case)]

extern crate merfishtools;
extern crate argparse;


use std::io::{stdout, stderr};
use std::str::FromStr;
use argparse::{ArgumentParser, Store, List, StoreTrue};

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

pub mod model;
pub mod io;
//pub mod cli;
/*
#[derive(Debug)]
enum Command {
    Exp,
    Diffexp,
    None
}


impl FromStr for Command {
    type Err = ();

    fn from_str(src: &str) -> Result<Command, ()> {
        return match src {
            "exp"     => Ok(Command::Exp),
            "diffexp" => Ok(Command::Diffexp),
            _            => Err(()),
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
        Command::Exp     => exp(args),
        Command::Diffexp => diffexp(args),
        Command::None       => {
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
    let mut threads = 1;

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"For given MERFISH data, calculate expressions for each feature (e.g. gene) in each cell.
Results are provided as PMF (probability mass function) in columns "Experiment, Cell, Feature, Expression, Posterior Probability".
Example: "merfishtools exp < data.txt > expression.txt""#
        );

        ap.refer(&mut N).add_option(&["-N"], Store, "Number of bits in readout (i.e. number of hybridization rounds).");
        ap.refer(&mut m).add_option(&["-m"], Store, "Number of 1-bits in readout.");
        ap.refer(&mut p0).add_option(&["--p0"], Store, "Probability of 0->1 error.");
        ap.refer(&mut p1).add_option(&["--p1"], Store, "Probability of 1->0 error.");
        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        parse_args_or_exit(&ap, args);
    }
    cli::expression(N, m, p0, p1, threads);
}


fn diffexp(args: Vec<String>) {
    let mut threads = 1;
    let mut group1_path = "".to_string();
    let mut group2_path = "".to_string();

    {
        let mut ap = ArgumentParser::new();
        ap.set_description(
r#"For given MERFISH data, calculate differentially expressed features (e.g. genes) between groups of cells given as separate input data.
Results are provided as PMF (probability mass function) in columns "Feature, Foldchange, Posterior Probability".
Example: "merfishtools diffexp data1.txt data2.txt > diffexp.txt""#
        );

        ap.refer(&mut threads)
          .add_option(&["--threads", "-t"], Store, "Number of threads to use.");
        ap.refer(&mut group1_path).required()
          .add_argument("group1", Store, "Group of cells.");
        ap.refer(&mut group2_path).required()
          .add_argument("group2", Store, "Group of cells.");
        parse_args_or_exit(&ap, args);
    }
    cli::differential_expression(&group1_path, &group2_path, threads);
}


fn parse_args_or_exit(ap: &ArgumentParser, args: Vec<String>) {
    match ap.parse(args, &mut stdout(), &mut stderr()) {
        Ok(()) =>  {}
        Err(x) => {
            std::process::exit(x);
        }
    }
}
*/
