#![allow(non_snake_case)]

use std;

use itertools::Itertools;
use simple_parallel;
use crossbeam;
use csv;

use bio::stats::logprobs::Prob;

use io;
use model;


pub fn expression(N: u8, m: u8, p0: Prob, p1: Prob, threads: usize) {
    let readout_model = model::Readout::new(N, m, p0, p1);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());
    let mut writer = io::pmf::expression::Writer::from_writer(std::io::stdout());

    let records = reader.records().map(
        |res| res.ok().expect("Error reading record.")
    ).group_by(|rec| {
        (rec.experiment, rec.cell_id, rec.feature.clone())
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
            writer.write(experiment, cell, &feature, &pmf);
        }
    });
}


pub fn differential_expression(group1_path: &str, group2_path: &str, threads: usize) {
    let mut reader1 = io::pmf::expression::Reader::from_file(group1_path);
    let mut reader2 = io::pmf::expression::Reader::from_file(group2_path);
    let mut writer = io::pmf::foldchange::Writer::from_writer(std::io::stdout());

    let group1 = reader1.pmfs();
    let group2 = reader2.pmfs();

    let features = group1.features().filter(|f| group2.contains_feature(f)).cloned().collect_vec();

    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for (feature, pmf) in pool.map(scope, features, |feature| {
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

            writer.write(&feature, &pmf);
        }
    });
}
