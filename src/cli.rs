#![allow(non_snake_case)]

use std;
use std::collections;

use itertools::Itertools;
use simple_parallel;
use crossbeam;
use regex::Regex;
use csv;

use bio::stats::logprobs::Prob;

use io;
use model;
use model::foldchange::LogFC;


pub struct Selection {
    pub expmnt: String,
    pub cell: String
}


#[derive(Debug)]
struct Counts {
    exact: u32,
    corrected: u32
}


pub fn stats(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8, codewords: u32, neighbors: u32) {
    let readout_model = model::Readout::new(N, m, p0, p1, dist, codewords, neighbors);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());

    let mut total_counts = collections::HashMap::new();
    for res in reader.records() {
        let record = res.unwrap();
        let total_count = total_counts.entry(record.cell_id).or_insert(0);
        *total_count += 1;
    }

    let mut writer = csv::Writer::from_writer(std::io::stdout()).delimiter(b'\t');
    writer.write(["cell", "total_expr_ev"].iter()).unwrap();

    for (cell, total_count) in total_counts {
        writer.write([cell, format!("{}", readout_model.expected_total(total_count))].iter()).unwrap();
    }
}


pub fn expression(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8, codewords: u32, neighbors: u32, estimate_path: Option<String>, threads: usize, cells: &str) {
    let readout_model = model::Readout::new(N, m, p0, p1, dist, codewords, neighbors);
    let mut reader = io::merfishdata::Reader::from_reader(std::io::stdin());
    let mut pmf_writer = io::pmf::expression::Writer::from_writer(std::io::stdout());
    let mut est_writer = estimate_path.map(|path| io::estimation::expression::Writer::from_file(path));

    let cells = Regex::new(cells).ok().expect("Invalid regular expression for cells.");

    let mut counts = collections::HashMap::new();
    let mut features = collections::HashSet::new();
    let mut count_total = 0;
    for record in reader.records().filter_map(|res| {
            let rec = res.unwrap();
            features.insert(rec.feature.clone());
            if cells.is_match(&rec.cell_id) {
                Some(rec)
            }
            else {
                None
            }
        }
    ) {
        count_total += 1;
        let cell_counts = counts.entry(record.cell_id).or_insert_with(collections::HashMap::new);
        let feature_counts = cell_counts.entry(record.feature).or_insert(Counts{ exact: 0, corrected: 0});
        if record.hamming_dist == 0 {
            feature_counts.exact += 1;
        }
        else if record.hamming_dist == 1 {
            feature_counts.corrected += 1;
        }
        else {
            panic!("Hamming distance of greater than 1 is unsupported at the moment.")
        }
    }
    for (_, cell_counts) in counts.iter_mut() {
        for feature in features.iter() {
            cell_counts.entry(feature.clone()).or_insert(Counts{ exact: 0, corrected: 0});
        }
    }

    let mut pool = simple_parallel::Pool::new(threads);
    crossbeam::scope(|scope| {
        for (_, (cell, pmfs)) in pool.unordered_map(scope, counts.into_iter(), |(cell, counts)| {
            let pmfs = counts.into_iter().map(|(feature, count)| {
                let pmf = model::expression::pmf(count.exact + count.corrected, count.exact, count_total, &readout_model);
                /*if feature == "COL5A1" {
                    debug!("{:?} {:?}", count, pmf);
                }*/
                (feature, pmf)
            }).collect_vec();
            (cell, pmfs)
        }) {
            for (feature, pmf) in pmfs {
                pmf_writer.write(&cell, &feature, &pmf);

                if let Some(ref mut est_writer) = est_writer {
                    est_writer.write(
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
