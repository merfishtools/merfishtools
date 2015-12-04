#![allow(non_snake_case)]

use std;
use std::collections;

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


struct Counts {
    exact: u32,
    corrected: u32
}


pub fn expression(N: u8, m: u8, p0: Prob, p1: Prob, dropout_rate: Prob, estimate_path: Option<String>, threads: usize, experiments: &str, cells: &str) {
    let readout_model = model::Readout::new(N, m, p0, p1, dropout_rate);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());
    let mut pmf_writer = io::pmf::expression::Writer::from_writer(std::io::stdout());
    let mut est_writer = estimate_path.map(|path| io::estimation::expression::Writer::from_file(path));

    let experiments = Regex::new(experiments).ok().expect("Invalid regular expression for experiments.");
    let cells = Regex::new(cells).ok().expect("Invalid regular expression for cells.");

    let mut counts = collections::HashMap::new();
    let mut features = collections::HashSet::new();
    for record in reader.records().filter_map(|res| {
            let rec = res.unwrap();
            features.insert(rec.feature.clone());
            if experiments.is_match(&rec.experiment) && cells.is_match(&rec.cell_id) {
                Some(rec)
            }
            else {
                None
            }
        }
    ) {
        let cell_counts = counts.entry((record.experiment, record.cell_id)).or_insert_with(collections::HashMap::new);
        let feature_counts = cell_counts.entry(record.feature).or_insert(Counts{ exact: 0, corrected: 0});
        if record.exact_match == 1 {
            feature_counts.exact += 1;
        }
        else {
            feature_counts.corrected += 1;
        }
    }
    for (_, cell_counts) in counts.iter_mut() {
        for feature in features.iter() {
            cell_counts.entry(feature.clone()).or_insert(Counts{ exact: 0, corrected: 0});
        }
    }


    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for (_, (experiment, cell, pmfs)) in pool.unordered_map(scope, counts.into_iter(), |((experiment, cell), counts)| {
            let pmfs = counts.into_iter().map(|(feature, count)| {
                (feature, model::expression::pmf(count.exact + count.corrected, count.corrected, &readout_model))
            }).collect_vec();
            (experiment, cell, pmfs)
        }) {
            for (feature, pmf) in pmfs {
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
        for (_, (feature, pmf)) in pool.unordered_map(scope, features, |feature| {
            info!("Calculating {}.", feature);
            let g1 = group1.get(&feature).unwrap();
            let g2 = group2.get(&feature).unwrap();
            let pmf = model::foldchange::pmf(
                &model::expressionset::pmf(&g1),
                &model::expressionset::pmf(&g2)
            );
            (
                feature,
                pmf
            )
        }) {
            if pmf.iter().any(|fc| fc.prob.is_nan()) {
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
