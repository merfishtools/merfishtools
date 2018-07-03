// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(non_snake_case)]

use std;
use std::collections;
use std::io::prelude::*;
use std::error::Error;

use csv;
use itertools::Itertools;
use cue;
use regex::Regex;

use bio::stats::{Prob, LogProb};

use io;
use model;
use model::foldchange::LogFC;
use model::cv::CV;
use codebook;
use model::readout::Counts;
use error_rates;


pub struct Selection {
    pub expmnt: String,
    pub cell: String
}


/*pub fn stats(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8, codebook_path: &str) {
    let codebook = io::codebook::Reader::from_file(codebook_path, dist).unwrap().codebook();
    let model = model::readout::new_model(N, m, p0, p1, dist, codebook);
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
        writer.write([cell, format!("{}", model.expected_total(total_count))].iter()).unwrap();
    }
}*/


/// Estimate expressions.
pub fn expression(
    p0: Vec<Prob>,
    p1: Vec<Prob>,
    codebook_path: &str,
    estimate_path: Option<&str>,
    stats_path: Option<&str>,
    threads: usize,
    cells: &str,
    window_width: u32,
    seed: usize
) {
    let codebook = io::codebook::Codebook::from_file(codebook_path).unwrap();
    let mut reader = io::merfishdata::tsv::Reader::from_reader(std::io::stdin());
    let mut cdf_writer = io::cdf::expression::Writer::from_writer(std::io::stdout());
    let mut est_writer = estimate_path.map(|path| io::estimation::expression::Writer::from_file(path));

    let mut stats_writer = stats_path.map(|path| csv::Writer::from_file(path).unwrap().delimiter(b'\t'));
    if let Some(ref mut writer) = stats_writer {
        writer.write(["cell", "noise-rate"].iter()).unwrap();
    }

    let cells = Regex::new(cells).ok().expect("Invalid regular expression for cells.");

    let mut counts = collections::HashMap::new();
    for record in reader.records().filter_map(|res| {
            let rec = res.unwrap();
            // consider record if it is contained in the codebook and has a valid cell id
            if cells.is_match(&rec.cell_id) && codebook.contains(&rec.feature) {
                Some(rec)
            }
            else {
                None
            }
        }
    ) {
        let cell_counts = counts.entry(record.cell_id).or_insert_with(collections::HashMap::new);
        let feature_counts = cell_counts.entry(record.feature).or_insert(Counts{ exact: 0, mismatch: 0});
        if record.hamming_dist == 0 {
            feature_counts.exact += 1;
        }
        else if record.hamming_dist == 1 {
            feature_counts.mismatch += 1;
        }
        else {
            panic!("Hamming distance of greater than 1 is unsupported at the moment.")
        }
    }
    // add missing features from codebook
    for (_, cell_counts) in counts.iter_mut() {
        for feature in codebook.features() {
            cell_counts.entry(feature.clone()).or_insert(Counts{ exact: 0, mismatch: 0});
        }
    }

    cue::pipeline(
        "exp",
        threads,
        counts.into_iter(),
        |(cell, counts)| {
            // Generate joint model.
            let mut model = model::readout::JointModel::new(
                counts.iter(), &p0, &p1, &codebook, window_width, seed
            );
            model.expectation_maximization(&cell);

            // Calculate CDF for all features.
            let cdfs = counts.into_iter().filter_map(|(feature, _)| {
                let feature_id = codebook.get_id(&feature);
                if codebook.record(feature_id).expressed() {
                    let (cdf, map_estimate) = model::expression::cdf(
                        feature_id, &mut model
                    );

                    Some((feature, cdf, map_estimate))
                } else { None }
            }).collect_vec();
            (cell, cdfs, model.noise_rate())
        },
        |(cell, cdfs, noise_rate)| {
            for (feature, cdf, map_estimate) in cdfs {
                cdf_writer.write(&cell, &feature, &cdf);

                if let Some(ref mut est_writer) = est_writer {
                    est_writer.write(
                        &cell,
                        &feature,
                        map_estimate,
                        cdf.credible_interval(0.95).expect("bug: empty CDF")
                    );
                }
            }
            if let Some(ref mut stats_writer) = stats_writer {
                stats_writer.write([&cell, &format!("{:.4}", noise_rate)].iter()).unwrap();
            }
        }
    );
}


/// Estimate differential expression over two conditions via fold changes.
pub fn differential_expression(group1_path: &str, group2_path: &str, pmf_path: Option<&str>, max_fc: LogFC, pseudocounts: f64, threads: usize) {
    assert!(pseudocounts > 0.0, "Pseudocounts must be > 0.0 for calculating fold changes.");

    let mut reader1 = io::cdf::expression::Reader::from_file(group1_path).expect("Invalid input for group 1.");
    let mut reader2 = io::cdf::expression::Reader::from_file(group2_path).expect("Invalid input for group 2.");
    let mut cdf_writer = pmf_path.map(|path| io::cdf::diffexp::Writer::from_file(path, "log2fc"));
    let mut est_writer = io::estimation::differential_expression::Writer::from_writer(std::io::stdout(), "log2fc");

    let group1 = reader1.cdfs();
    let group2 = reader2.cdfs();

    let features = group1.features().filter(|f| group2.contains_feature(f)).cloned().collect_vec();

    let mut estimates = Vec::new();

    cue::pipeline(
        "diffexp",
        threads,
        features.iter(),
        |feature| {
            info!("Calculating {}.", feature);
            let g1 = group1.get(&feature).unwrap();
            let g2 = group2.get(&feature).unwrap();
            let cdf = model::foldchange::cdf(
                &model::expressionset::cdf(&g1, pseudocounts),
                &model::expressionset::cdf(&g2, pseudocounts)
            );
            (feature, cdf)
        },
        |(feature, cdf)| {
            if cdf.iter().any(|e| e.prob.is_nan()) {
                warn!("A CDF probability for feature {} cannot be estimated because counts are too high. It will be reported as NaN.", feature);
            }

            if let Some(ref mut cdf_writer) = cdf_writer {
                cdf_writer.write(&feature, &cdf);
            }
            estimates.push((feature, model::diffexp::estimate(&cdf, max_fc)));
        }
    );

    // calculate FDR and write estimates to STDOUT
    estimates.sort_by(|&(_, ref a), &(_, ref b)| {
        a.differential_expression_pep.partial_cmp(&b.differential_expression_pep).unwrap()
    });
    let expected_fds = LogProb::ln_cumsum_exp(
        estimates.iter().map(|&(_, ref a)| a.differential_expression_pep)
    ).collect_vec();
    for (i, ((feature, estimate), fd)) in estimates.into_iter().zip(expected_fds).enumerate() {
        let fdr = LogProb(*fd - (i as f64 + 1.0).ln());
        est_writer.write(
            &feature,
            estimate.differential_expression_pep,
            fdr,
            estimate.differential_expression_bf,
            estimate.map,
            estimate.credible_interval
        );
    }
}


/// Estimate differential expression over multiple conditions via the coefficient of variation.
pub fn multi_differential_expression(group_paths: &[&str], pmf_path: Option<&str>, max_cv: CV, pseudocounts: f64, threads: usize) {
    let groups = group_paths.iter().enumerate().map(|(i, path)| {
        io::cdf::expression::Reader::from_file(path).expect(&format!("Invalid input for group {}.", i))
                                                    .cdfs()
    }).collect_vec();
    let mut cdf_writer = pmf_path.map(|path| io::cdf::diffexp::Writer::from_file(path, "cv"));
    let mut est_writer = io::estimation::differential_expression::Writer::from_writer(std::io::stdout(), "cv");

    // TODO take feature intersection or warn if features are not the same
    let mut features = Vec::new();
    for feature in groups[0].features() {
        let mut common = true;
        for group in groups[1..].iter() {
            common &= group.contains_feature(feature);
        }
        if common {
            features.push(feature.clone());
        }
    }

    let mut estimates = Vec::new();

    cue::pipeline(
        "multidiffexp",
        threads,
        features.iter(),
        |feature| {
            info!("Calculating {}.", feature);
            let cdfs = groups.iter().map(|group| {
                model::expressionset::cdf(
                    group.get(&feature).expect("Missing feature."),
                    pseudocounts
                )
            }).collect_vec();
            let cdf = model::cv::cdf(&cdfs);
            (feature, cdf)
        },
        |(feature, cdf)| {
            if cdf.iter().any(|e| e.prob.is_nan()) {
                warn!("A CDF probability for feature {} cannot be estimated due to numerical problems. It will be reported as NaN.", feature);
            }

            if let Some(ref mut cdf_writer) = cdf_writer {
                cdf_writer.write(&feature, &cdf);
            }
            estimates.push((feature, model::diffexp::estimate(&cdf, max_cv)));
        }
    );

    // calculate FDR and write estimates to STDOUT
    estimates.sort_by(|&(_, ref a), &(_, ref b)| {
        a.differential_expression_pep.partial_cmp(&b.differential_expression_pep).unwrap()
    });
    let expected_fds = LogProb::ln_cumsum_exp(
        estimates.iter().map(|&(_, ref a)| a.differential_expression_pep)
    ).collect_vec();
    for (i, ((feature, estimate), fd)) in estimates.into_iter().zip(expected_fds).enumerate() {
        let fdr = LogProb(*fd - (i as f64 + 1.0).ln());
        est_writer.write(
            &feature,
            estimate.differential_expression_pep,
            fdr,
            estimate.differential_expression_bf,
            estimate.map,
            estimate.credible_interval
        );
    }
}


pub fn gen_codebook(
    words: &[codebook::Word],
    not_expressed_pattern: Option<&str>
) -> Result<(), Box<Error>> {
    let not_expressed_re = if not_expressed_pattern.is_some() {
        Some(Regex::new(not_expressed_pattern.unwrap())?)
    } else {
        None
    };

    let stdin = std::io::stdin();
    let mut reader = stdin.lock().lines();
    let mut writer = csv::Writer::from_writer(std::io::stdout()).delimiter(b'\t');
    let mut words = words.iter();
    // TODO add expressed column and handle misidentification probes
    writer.write(["feat", "codeword", "expressed"].iter()).unwrap();

    for i in 1.. {
        match (reader.next(), words.next()) {
            (Some(feature), Some(w)) => {
                let feature = feature?;
                if feature.len() == 0 {
                    // TODO proper error handling
                    panic!("Empty feature found. All features provided at STDIN have to be non-empty.");
                }

                let expressed = !not_expressed_re.as_ref().map_or_else(
                    || false, |re| re.is_match(&feature)
                );

                writer.write([
                    feature,
                    format!("{:?}", w),
                    format!("{}", if expressed { "1" } else { "0" })
                ].into_iter())?;
            },
            (None, Some(_)) => break,
            (Some(_), None) => {
                error!("Not enough code words. Generated codebook for the first {} transcripts.", i);
                break
            },
            (None, None) => break
        }
    }
    Ok(())
}


pub fn estimate_error_rates(
    codebook: &str
) -> Result<(), Box<Error>> {

    let codebook = io::codebook::Codebook::from_file(codebook).unwrap();
    let mut readouts = csv::Reader::from_reader(std::io::stdin()).delimiter(b'\t');

    let (p0, p1) = error_rates::estimate(
        &codebook,
        readouts.decode().map(|rec| {
            let (cell, feat, readout): (String, String, String) = rec.unwrap();
            (feat, io::codebook::parse_codeword(readout.as_bytes()))
        })
    );

    let mut writer = csv::Writer::from_writer(std::io::stdout()).delimiter(b'\t');
    writer.write(["pos", "p0", "p1"].iter())?;
    for (i, (p0, p1)) in p0.into_iter().zip(p1.into_iter()).enumerate() {
        writer.write([format!("{}", i), format!("{}", *p0), format!("{}", *p1)].iter())?;
    }
    Ok(())
}
