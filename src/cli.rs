use std;

use itertools::Itertools;
use simple_parallel;
use crossbeam;

use bio::stats::logprobs::Prob;

use io;
use model;

pub fn expression(N: u8, m: u8, p0: Prob, p1: Prob, threads: usize) {
    let readout_model = model::Readout::new(N, m, p0, p1);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());
    let mut writer = io::pmf::Writer::from_writer(std::io::stdout());

    let mut records = reader.records().map(
        |res| res.ok().expect("Error reading record.")
    ).group_by(|rec| {
        (rec.experiment, rec.cell_id, rec.feature.clone())
    });

    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for ((experiment, cell, feature), expression) in pool.map(scope, records, |((experiment, cell, feature), readouts)| {
            let count = readouts.len();
            let count_exact = readouts.iter().filter(|r| r.exact_match == 1).count();
            ((experiment, cell, feature), model::Expression::new(count as u32, count_exact as u32, &readout_model))
        }) {
            let mut record = io::pmf::Record{
                feature: io::expression::Record{
                    experiment: experiment,
                    cell: cell,
                    feature: feature.clone()
                },
                value: 0,
                prob: 0.0
            };

            for (value, prob) in expression.pmf() {
                record.value = value;
                record.prob = prob;
                writer.write(&record);
            }
        }
    });
}


pub fn differential_expression(group1_path: &str, group2_path: &str, threads: usize) {
    let reader1 = io::pmf::Reader::from_file(group1_path);
    let reader2 = io::pmf::Reader::from_file(group2_path);
}
