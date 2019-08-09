// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

#![allow(non_snake_case)]

use std;
use std::collections;
use std::io::prelude::*;

use bio::stats::{LogProb, Prob};
use csv;
use cue;
use failure::Error;
use itertools::Itertools;
use regex::Regex;

use crate::codebook;
use crate::error_rates;
use crate::io;
use crate::io::merfishdata::{self, MerfishRecord, Reader};
use crate::model::bayes;
use crate::model::bayes::cv::CV;
use crate::model::bayes::foldchange::LogFC;
use crate::model::bayes::readout::Counts;

pub struct Selection {
    pub expmnt: String,
    pub cell: String,
}

/*pub fn stats(N: u8, m: u8, p0: Prob, p1: Prob, dist: u8, codebook_path: &str) {
    let codebook = io::codebook::Reader::from_file(codebook_path, dist).unwrap().codebook();
    let model = bayes::readout::new_model(N, m, p0, p1, dist, codebook);
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

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct ExpressionJ {
    p0: Vec<Prob>,
    p1: Vec<Prob>,
    codebook_path: String,
    estimate_path: Option<String>,
    stats_path: Option<String>,
    threads: usize,
    cells: Regex,
    window_width: u32,
    seed: u64,
    #[builder(setter(skip))]
    counts: collections::HashMap<String, collections::HashMap<String, Counts>>,
}

pub trait Expression {
    fn codebook_path(&self) -> &str;
    fn cells(&self) -> &Regex;
    fn infer(&mut self) -> Result<(), Error>;
}

impl Expression for ExpressionJ {
    fn codebook_path(&self) -> &str {
        &self.codebook_path
    }

    fn cells(&self) -> &Regex {
        &self.cells
    }

    /// Estimate expressions.
    fn infer(&mut self) -> Result<(), Error> {
        let codebook = io::codebook::Codebook::from_file(&self.codebook_path()).unwrap();
        let mut cdf_writer = io::cdf::expression::Writer::from_writer(std::io::stdout());
        let mut est_writer = self
            .estimate_path
            .as_ref()
            .map(io::estimation::expression::Writer::from_file);

        let mut stats_writer = self.stats_path.as_ref().map(|path| {
            csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(path)
                .unwrap()
        });
        if let Some(ref mut writer) = stats_writer {
            writer.serialize(["cell", "noise-rate"]).unwrap();
        }

        cue::pipeline(
            "exp",
            self.threads,
            self.counts.iter(),
            |(cell, counts)| {
                // Generate joint model.
                let mut model = bayes::readout::JointModel::new(
                    counts.iter(),
                    &self.p0,
                    &self.p1,
                    &codebook,
                    self.window_width,
                    self.seed,
                );
                model.expectation_maximization(&cell);

                // Calculate CDF for all features.
                let cdfs = counts
                    .into_iter()
                    .filter_map(|(feature, _)| {
                        let feature_id = codebook.get_id(&feature);
                        if codebook.record(feature_id).expressed() {
                            let (cdf, map_estimate) =
                                bayes::expression::cdf(feature_id, &mut model);

                            Some((feature, cdf, map_estimate))
                        } else {
                            None
                        }
                    })
                    .collect_vec();
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
                            &cdf.credible_interval(0.95).expect("bug: empty CDF"),
                        );
                    }
                }
                if let Some(ref mut stats_writer) = stats_writer {
                    stats_writer
                        .serialize([&cell, &format!("{:.4}", noise_rate)])
                        .unwrap();
                }
            },
        );

        Ok(())
    }
}

impl ExpressionJ {
    pub fn load_counts<'a, R>(
        &mut self,
        reader: &'a mut R,
        _format: crate::io::merfishdata::Format,
    ) -> Result<(), Error>
    where
        R: io::merfishdata::Reader<'a>,
    {
        let codebook = io::codebook::Codebook::from_file(&self.codebook_path())?;

        let mut counts = collections::HashMap::new();
        for res in reader.records() {
            let record = res?;
            if !self.cells().is_match(&record.cell_name())
                && codebook.contains(&record.feature_name())
            {
                continue;
            }

            let cell_counts = counts
                .entry(record.cell_name())
                .or_insert_with(collections::HashMap::new);
            let feature_counts = cell_counts.entry(record.feature_name()).or_insert(Counts {
                exact: 0,
                corrected: 0,
                uncorrected: 0,
            });
            if record.hamming_dist() == 0 {
                feature_counts.exact += 1;
            } else if record.hamming_dist() == 1 {
                feature_counts.corrected += 1;
            } else {
                feature_counts.uncorrected += 1;
                panic!("Hamming distance of greater than 1 is unsupported at the moment.")
            }
        }
        // add missing features from codebook
        for cell_counts in counts.values_mut() {
            for feature in codebook.features() {
                cell_counts.entry(feature.clone()).or_insert(Counts {
                    exact: 0,
                    corrected: 0,
                    uncorrected: 0,
                });
            }
        }
        self.counts = counts;
        Ok(())
    }
}

/// Estimate differential expression over two conditions via fold changes.
pub fn differential_expression(
    group1_path: &str,
    group2_path: &str,
    pmf_path: Option<&str>,
    max_fc: LogFC,
    pseudocounts: f64,
    threads: usize,
) -> Result<(), Error> {
    assert!(
        pseudocounts > 0.0,
        "Pseudocounts must be > 0.0 for calculating fold changes."
    );

    let mut reader1 = io::cdf::expression::Reader::from_file(group1_path)?;
    let mut reader2 = io::cdf::expression::Reader::from_file(group2_path)?;
    let mut cdf_writer = pmf_path.map(|path| io::cdf::diffexp::Writer::from_file(path, "log2fc"));
    let mut est_writer =
        io::estimation::differential_expression::Writer::from_writer(std::io::stdout(), "log2fc");

    let group1 = reader1.cdfs();
    let group2 = reader2.cdfs();

    let features = group1
        .features()
        .filter(|f| group2.contains_feature(f))
        .cloned()
        .collect_vec();

    let mut estimates = Vec::new();

    cue::pipeline(
        "diffexp",
        threads,
        features.iter(),
        |feature| {
            info!("Calculating {}.", feature);
            let g1 = group1.get(&feature).unwrap();
            let g2 = group2.get(&feature).unwrap();
            let cdf = bayes::foldchange::cdf(
                &bayes::expressionset::cdf(&g1, pseudocounts),
                &bayes::expressionset::cdf(&g2, pseudocounts),
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
            estimates.push((feature, bayes::diffexp::estimate(&cdf, max_fc)));
        },
    );

    // calculate FDR and write estimates to STDOUT
    estimates.sort_by(|&(_, ref a), &(_, ref b)| {
        a.differential_expression_pep
            .partial_cmp(&b.differential_expression_pep)
            .unwrap()
    });
    let expected_fds = LogProb::ln_cumsum_exp(
        estimates
            .iter()
            .map(|&(_, ref a)| a.differential_expression_pep),
    )
    .collect_vec();

    for (i, ((feature, estimate), fd)) in estimates.into_iter().zip(expected_fds).enumerate() {
        let fdr = LogProb(*fd - (i as f64 + 1.0).ln());
        est_writer.write(
            &feature,
            estimate.differential_expression_pep,
            fdr,
            estimate.differential_expression_bf,
            estimate.map,
            estimate.credible_interval,
        );
    }

    Ok(())
}

/// Estimate differential expression over multiple conditions via the coefficient of variation.
pub fn multi_differential_expression(
    group_paths: &[&str],
    pmf_path: Option<&str>,
    max_cv: CV,
    pseudocounts: f64,
    threads: usize,
) -> Result<(), Error> {
    let mut groups = Vec::new();
    for path in group_paths {
        let mut reader = io::cdf::expression::Reader::from_file(path)?;
        groups.push(reader.cdfs());
    }
    let mut cdf_writer = pmf_path.map(|path| io::cdf::diffexp::Writer::from_file(path, "cv"));
    let mut est_writer =
        io::estimation::differential_expression::Writer::from_writer(std::io::stdout(), "cv");

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
            let cdfs = groups
                .iter()
                .map(|group| {
                    bayes::expressionset::cdf(
                        group.get(&feature).expect("Missing feature."),
                        pseudocounts,
                    )
                })
                .collect_vec();
            let cdf = bayes::cv::cdf(&cdfs);
            (feature, cdf)
        },
        |(feature, cdf)| {
            if cdf.iter().any(|e| e.prob.is_nan()) {
                warn!("A CDF probability for feature {} cannot be estimated due to numerical problems. It will be reported as NaN.", feature);
            }

            if let Some(ref mut cdf_writer) = cdf_writer {
                cdf_writer.write(&feature, &cdf);
            }
            estimates.push((feature, bayes::diffexp::estimate(&cdf, max_cv)));
        },
    );

    // calculate FDR and write estimates to STDOUT
    estimates.sort_by(|&(_, ref a), &(_, ref b)| {
        a.differential_expression_pep
            .partial_cmp(&b.differential_expression_pep)
            .unwrap()
    });
    let expected_fds = LogProb::ln_cumsum_exp(
        estimates
            .iter()
            .map(|&(_, ref a)| a.differential_expression_pep),
    )
    .collect_vec();
    for (i, ((feature, estimate), fd)) in estimates.into_iter().zip(expected_fds).enumerate() {
        let fdr = LogProb(*fd - (i as f64 + 1.0).ln());
        est_writer.write(
            &feature,
            estimate.differential_expression_pep,
            fdr,
            estimate.differential_expression_bf,
            estimate.map,
            estimate.credible_interval,
        );
    }

    Ok(())
}

pub fn gen_codebook(
    words: &[codebook::Word],
    not_expressed_pattern: Option<&str>,
) -> Result<(), Error> {
    let not_expressed_re = if not_expressed_pattern.is_some() {
        Some(Regex::new(not_expressed_pattern.unwrap())?)
    } else {
        None
    };

    let stdin = std::io::stdin();
    let mut reader = stdin.lock().lines();
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(std::io::stdout());
    let mut words = words.iter();
    // TODO add expressed column and handle misidentification probes
    writer.serialize(["feat", "codeword", "expressed"]).unwrap();

    for i in 1.. {
        match (reader.next(), words.next()) {
            (Some(feature), Some(w)) => {
                let feature = feature?;
                if feature.is_empty() {
                    // TODO proper error handling
                    panic!(
                        "Empty feature found. All features provided at STDIN have to be non-empty."
                    );
                }

                let expressed = !not_expressed_re
                    .as_ref()
                    .map_or_else(|| false, |re| re.is_match(&feature));

                writer.serialize([
                    feature,
                    format!("{:?}", w),
                    if expressed { "1" } else { "0" }.to_string(),
                ])?;
            }
            (None, Some(_)) => break,
            (Some(_), None) => {
                error!(
                    "Not enough code words. Generated codebook for the first {} transcripts.",
                    i
                );
                break;
            }
            (None, None) => break,
        }
    }
    Ok(())
}

pub fn estimate_error_rates(raw_data: &str, codebook: &str) -> Result<(), Error> {
    let codebook = io::codebook::Codebook::from_file(codebook).unwrap();

    let (p0, p1) = if let merfishdata::Format::Binary = merfishdata::Format::from_path(raw_data) {
        let mut reader = merfishdata::binary::BinaryReader::from_file(raw_data)?;
        error_rates::estimate(
            &codebook,
            reader.records().map(|rec| {
                let rec = rec.unwrap();
                (codebook.get_name(rec.feature_id()).unwrap(), rec.readout_bitvec())
            }),
        )
    } else {
        let mut readouts = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(raw_data)?;

        error_rates::estimate(
            &codebook,
            readouts.deserialize().map(|rec| {
                let (_cell, feat, readout): (String, String, String) = rec.unwrap();
                (feat, io::codebook::parse_codeword(readout.as_bytes()))
            }),
        )
    };

    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(std::io::stdout());
    writer.write_record(&["pos", "p0", "p1"])?;
    for (i, (p0, p1)) in p0.into_iter().zip(p1.into_iter()).enumerate() {
        writer.write_record(&[format!("{}", i), format!("{}", *p0), format!("{}", *p1)])?;
    }
    Ok(())
}
