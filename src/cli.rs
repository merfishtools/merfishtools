#![allow(non_snake_case)]

use std;

use itertools::Itertools;
use simple_parallel;
use crossbeam;
use regex::Regex;

use bio::stats::logprobs::Prob;

use io;
use model;
use model::foldchange::LogFC;


pub struct Selection {
    pub expmnt: String,
    pub cell: String
}


pub fn expression(N: u8, m: u8, p0: Prob, p1: Prob, estimate_path: Option<String>, threads: usize, experiments: &str, cells: &str) {
    let readout_model = model::Readout::new(N, m, p0, p1);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());
    let mut pmf_writer = io::pmf::expression::Writer::from_writer(std::io::stdout());
    let mut est_writer = estimate_path.map(|path| io::estimation::expression::Writer::from_file(path));

    let experiments = Regex::new(experiments).ok().expect("Invalid regular expression for experiments.");
    let cells = Regex::new(cells).ok().expect("Invalid regular expression for cells.");

    let records = reader.records().filter_map(|res| {
            let rec = res.unwrap();
            if experiments.is_match(&rec.experiment) && cells.is_match(&rec.cell_id) {
                Some(rec)
            }
            else {
                None
            }
        }
    ).group_by(|rec| {
        (rec.experiment.clone(), rec.cell_id.clone(), rec.feature.clone())
    });

    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for ((experiment, cell, feature), pmf) in pool.map(scope, records, |((experiment, cell, feature), readouts)| {
            let count = readouts.len();
            let count_exact = readouts.iter().filter(|r| r.exact_match == 1).count();
            (
                (experiment, cell, feature),
                model::expression::pmf(count as u32, count_exact as u32, &readout_model)
            )
        }) {
            pmf_writer.write(&experiment, &cell, &feature, &pmf);

            if let Some(ref mut est_writer) = est_writer {
                est_writer.write(
                    &experiment,
                    &cell,
                    &feature,
                    pmf.expected_value(),
                    pmf.standard_deviation()
                );
            }
        }
    });
}


pub fn differential_expression(group1_path: &str, group2_path: &str, pmf_path: Option<String>, min_fc: LogFC, threads: usize) {
    let mut reader1 = io::pmf::expression::Reader::from_file(group1_path);
    let mut reader2 = io::pmf::expression::Reader::from_file(group2_path);
    let mut pmf_writer = pmf_path.map(|path| io::pmf::foldchange::Writer::from_file(path));
    let mut est_writer = io::estimation::differential_expression::Writer::from_writer(std::io::stdout());

    let group1 = reader1.pmfs();
    let group2 = reader2.pmfs();

    let features = group1.features().filter(|f| group2.contains_feature(f)).cloned().collect_vec();

    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for (feature, pmf) in pool.map(scope, features, |feature| {
            info!("Calculating {}.", feature);
            let pmf = model::foldchange::pmf(
                &model::expressionset::pmf(&group1.get(&feature).unwrap()),
                &model::expressionset::pmf(&group2.get(&feature).unwrap())
            );
            (
                feature,
                pmf
            )
        }) {
            if pmf.iter().any(|&(_, prob)| prob.is_nan()) {
                warn!("A PMF value for feature {} cannot be estimated because counts are too high. It will be reported as NaN.", feature);
            }

            est_writer.write(
                &feature,
                pmf.differential_expression_pep(min_fc),
                pmf.expected_value(),
                pmf.standard_deviation()
            );

            if let Some(ref mut pmf_writer) = pmf_writer {
                pmf_writer.write(&feature, &pmf);
            }
        }
    });
}
