use std::collections::HashMap;

use bio::stats::Prob;
use failure::Error;
use ndarray::prelude::*;
use rand;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use clap::arg_enum;

use crate::cli::Expression;
use crate::io::codebook::{Codebook, FeatureID};
use crate::io::merfishdata::MerfishRecord;
use crate::model::bayes::readout::Counts;
use crate::model::la::common::{Errors, Expr, ExprV, NUM_BITS, NUM_CODES};
use crate::model::la::expression::Mode::ExpressionThenErrors;
use crate::model::la::matrix::{CSR, csr_error_matrix, csr_successive_overrelaxation, error_dot, error_successive_overrelaxation};
use crate::model::la::problem::{objective, partial_objective};
use crate::io::merfishdata;
use crate::model::la::hamming::hamming_distance16;
use itertools::Itertools;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct ExpressionT {
    p0: Vec<Prob>,
    p1: Vec<Prob>,
    codebook_path: String,
    estimate_path: Option<String>,
    threads: usize,
    cells: Regex,
    bits: usize,
    max_hamming_distance: usize,
    mode: Mode,
    absolute: bool,
    seed: u64,
    #[builder(setter(skip))]
    corrected_counts: HashMap<String, HashMap<u16, usize>>,
    #[builder(setter(skip))]
    raw_counts: HashMap<String, HashMap<u16, usize>>,
}

impl Expression for ExpressionT {
    fn codebook_path(&self) -> &str {
        &self.codebook_path
    }

    fn cells(&self) -> &Regex {
        &self.cells
    }

    fn infer(&mut self) -> Result<(), Error> {
        let codebook = Codebook::from_file(&self.codebook_path()).unwrap();
        let codemap: HashMap<&String, usize> = codebook.features().map(|feature| {
            let feature_id = codebook.get_id(feature);
            let idx: usize = codebook.record(feature_id).codeword()
                .iter().rev().enumerate().map(|(i, bit)| (bit as usize) << i).sum();
            (feature, idx)
        }).collect();

        let raw_counts = &self.raw_counts;
        let corrected_counts = &self.corrected_counts;
        let max_hamming_distance = self.max_hamming_distance;
        let mode = &self.mode;

        let mut e = [[0.; NUM_BITS], [0.; NUM_BITS]];
        e[0].iter_mut().zip(self.p0.iter().cloned().map(f64::from))
            .for_each(|(v, p)| *v = p as f32);
        e[1].iter_mut().zip(self.p1.iter().cloned().map(f64::from))
            .for_each(|(v, p)| *v = p as f32);

        let mut x_est: Expr = Array1::<f32>::zeros(NUM_CODES);
        let mut y: Expr = Array1::<f32>::zeros(NUM_CODES);
        for (_cell, feature_counts) in corrected_counts {
            for (&barcode, &count) in feature_counts {
                x_est[barcode as usize] += count as f32;
            }
        }
        for (_cell, feature_counts) in raw_counts {
            for (&barcode, &count) in feature_counts {
                y[barcode as usize] += count as f32;
            }
        }
//        for (i, (&xx, &yy)) in x_est.iter().zip(y.iter()).enumerate() {
//            if xx != 0. || yy != 0. {
//                dbg!((i, (xx, yy)));
//            }
//        }
//        std::process::exit(0);
        let magnitude = y.sum();
        if !self.absolute {
            y /= magnitude;
            x_est /= magnitude;
        }
        let codewords = codemap.values().collect_vec();

        let yv = y.view();
        let mut x: Expr = x_est;

        match mode {
            Mode::ErrorsThenExpression => {
                for hamming_dist in 1..=max_hamming_distance {
                    println!("\n=== {:?} ===", hamming_dist);
                    println!("Estimating positional error probabilities via GD");
                    e = estimate_errors(x.view(), yv, &e, hamming_dist).expect("Failed estimating errors");
                    dbg!(&e);
                    println!("Constructing CSR transition matrix");
                    let mat = csr_error_matrix(&e, hamming_dist.max(2));
                    println!("Estimating true expression via SOR");
                    x = estimate_expression(&mat, x.view(), yv, hamming_dist, false).expect("Failed estimating expression");
                    x.iter_mut().filter(|&&mut v| v <= 0.).for_each(|v| *v = 0.);
                    if !self.absolute {
                        x /= x.sum();
                    }
                    dbg!(codewords.iter().map(|&i| y[*i]).collect::<Vec<f32>>());
                    dbg!(codewords.iter().map(|&i| x[*i]).collect::<Vec<f32>>());
                }
            }
            Mode::ExpressionThenErrors => {
                for hamming_dist in 1..=max_hamming_distance {
                    println!("\n=== {:?} ===", hamming_dist);
                    println!("Constructing CSR transition matrix");
                    let mat = csr_error_matrix(&e, hamming_dist.max(2));
                    println!("Estimating true expression via SOR");
                    x = estimate_expression(&mat, x.view(), yv, hamming_dist, false).expect("Failed estimating expression");
                    x.iter_mut().filter(|&&mut v| v <= 0.).for_each(|v| *v = 0.);
                    if !self.absolute {
                        x /= x.sum();
                    }
                    dbg!(codewords.iter().map(|&i| y[*i]).collect::<Vec<f32>>());
                    dbg!(codewords.iter().map(|&i| x[*i]).collect::<Vec<f32>>());
                    println!("Estimating positional error probabilities via GD");
                    e = estimate_errors(x.view(), yv, &e, hamming_dist).expect("Failed estimating errors");
                    dbg!(&e);
                }
            }
        };
        let mat = csr_error_matrix(&e, 4);
        self.raw_counts
            .par_iter()
            .for_each(|(cell, raw_countmap)| {
                let corrected_countmap = self.corrected_counts.get(cell).unwrap();
                let mut y: Expr = Expr::zeros((NUM_CODES, ));
                raw_countmap.iter().for_each(|(&barcode, &count)| {
                    y[barcode as usize] += count as f32;
                });
                let mut x: Expr = Expr::zeros((NUM_CODES, ));
                corrected_countmap.iter().for_each(|(&barcode, &count)| {
                    x[barcode as usize] += count as f32;
                });
                let magnitude_y = y.sum();
                let magnitude_x = x.sum();
                dbg!(magnitude_x);
                if !self.absolute {
                    y /= magnitude_y;
                    x /= magnitude_x;
                }
                x = estimate_expression(&mat, x.view(), y.view(), max_hamming_distance, false).expect("Failed estimating expression");
                x *= magnitude_x;
//                x.iter_mut().filter(|&&mut v| v <= 1.).for_each(|v| *v = 0.);
//                x.iter_mut().enumerate().filter(|(i, v)| corrected_countmap.contains_key(&(*i as u16))).for_each(|(_, v)| *v = 0.);
//                dbg!(codewords.iter().map(|&i| y[*i]).collect::<Vec<f32>>());
//                dbg!(codewords.iter().map(|&i| x[*i]).collect::<Vec<f32>>());
//                x.iter().enumerate().filter(|(i, &v)| v > 0.).for_each(|(i, v)| println!("{}\t{:016b}\t{}", cell, i, v));
                codewords.iter().map(|&i| (i, x[*i])).for_each(|(i, v)| println!("{}\t{:016b}\t{}", cell, i, v));
                println!();
            });
        Ok(())
    }
}

impl ExpressionT {
    pub fn load_counts<'a, R>(&mut self, reader: &'a mut R, format: merfishdata::Format) -> Result<(), Error>
        where R: crate::io::merfishdata::Reader<'a> {
        let codebook = &crate::io::codebook::Codebook::from_file(&self.codebook_path())?;

        let mut corrected_counts: HashMap<String, HashMap<u16, usize>> = HashMap::new();
        let mut raw_counts: HashMap<String, HashMap<u16, usize>> = HashMap::new();

        let mut rng = StdRng::seed_from_u64(self.seed);

        let codewords: Vec<u16> = codebook.records()
            .iter()
            .map(|r| r.codeword().iter().rev().enumerate().map(|(i, bit)| (bit as u16) << i).sum())
            .collect();

        for rec in reader.records() {
            let record = rec?;
            if !self.cells().is_match(&record.cell_name()) && codebook.contains(&record.feature_name()) {
                continue;
            }
            let barcode = record.barcode(Some(codebook));
            let mut raw_barcode = barcode;
            let mut uncorrected_barcode = barcode;
            match format {
                // simulated data is perturbed ground truth data,
                // i.e. no 1-bit error correction has been performed
                // ... so we do that now
                merfishdata::Format::Simulation => {
                    if record.hamming_dist() == 1 {
                        let closest: Vec<u16> = codewords.iter()
                            .map(|&cw| (cw, hamming_distance16(cw as usize, raw_barcode as usize)))
                            .filter(|(cw, d)| *d == 1)
                            .map(|(cw, _)| cw)
                            .collect();
                        if !closest.is_empty() {
                            raw_barcode = closest[0];
                        }
                    }
                }
                merfishdata::Format::TSV | merfishdata::Format::Binary => {
                    // if the readout has been corrected, reconstruct the/an un-corrected barcode
                    if record.hamming_dist() == 1 {
                        let error_bit = match record.error_bit() {
                            // the binary merfish format actually tells us at which position the error occurred
                            Some(bit) => bit,

                            // the tsv merfish format only states if a readout was bit-corrected or not
                            // we will hence un-correct the barcode at random (according to p0/p1 probabilities)
                            None => {
                                let weights: Vec<_> = (0..NUM_BITS)
                                    .map(|i| {
                                        let bit = (barcode >> i) & 1;
                                        let p = match bit {
                                            0 => &self.p0[i],
                                            1 => &self.p1[i],
                                            _ => panic!("bug: bit can only be 0 or 1")
                                        };
                                        p.0
                                    }).collect();
                                let wc = WeightedIndex::new(&weights).unwrap();
                                rng.sample(&wc) as u8
                            }
                        };
                        uncorrected_barcode ^= 1 << error_bit;
                    }
                }
            }

            // TODO build Map<Cell, Map<Barcode, Count>>
            let cell_counts_corrected = corrected_counts.entry(record.cell_name()).or_insert_with(HashMap::new);
            let feature_counts_corrected = cell_counts_corrected.entry(raw_barcode).or_insert(0);
            *feature_counts_corrected += 1;

            let cell_counts_raw = raw_counts.entry(record.cell_name()).or_insert_with(HashMap::new);
            let feature_counts_raw = cell_counts_raw.entry(uncorrected_barcode).or_insert(0);
            *feature_counts_raw += 1;
        }
        self.corrected_counts = corrected_counts;
        self.raw_counts = raw_counts;
        Ok(())
    }
}

fn _gradient_descent_hardcoded(e: &mut Errors,
                               y: ExprV,
                               x: ExprV,
                               target_error: f32,
                               alpha: f32,
                               max_hamming_distance: usize,
                               min_iter: usize,
                               max_iter: usize,
                               x_ind: &[usize]) -> (Errors, f32, usize) {
    let mut gradient: Errors = [[0.; NUM_BITS]; 2];
    let e0 = e;
    let alpha0 = alpha;
    let mult: f32 = 0.5;
    let c1: f32 = 1e-4;
    for i in 0..max_iter {
        println!("Calculating gradient...");
        (0..2).for_each(|kind| {
            (0..NUM_BITS).for_each(|pos| {
                let grad_value = partial_objective(x, y, &e0, pos, kind, max_hamming_distance, x_ind);
                gradient[kind][pos] = grad_value;
            });
        });
        println!("Estimating stepsize...");
        let o2 = objective(x, y, &e0, max_hamming_distance, x_ind);
        // gradv → descent direction: -partial_objective^T · partial_objective
        let gradv = -gradient.iter().flatten().map(|g: &f32| g.powi(2)).sum::<f32>();
        dbg!(&gradient);
        dbg!((o2, gradv));
        let mut alpha = alpha0;
        for _ in 0..17 {
            let mut e1 = *e0;
            for kind in 0..2 {
                for pos in 0..NUM_BITS {
                    let v = alpha * gradient[kind][pos];
                    e1[kind][pos] -= v;
                    if e1[kind][pos] < 0. {
                        e1[kind][pos] *= -1.;
                    }
                }
            }
            // e1 = e0 + alpha · gradv
            let o1 = objective(x, y, &e1, max_hamming_distance, x_ind);
            let o2 = o2 + c1 * alpha * gradv;
            // armijo rule →  f(e0 + alpha · gradv) <= f(e0) + c1 · alpha · gradv?
            if o1 <= o2 {
                break;
            }
            alpha *= mult;
        }
        for kind in 0..2 {
            for pos in 0..NUM_BITS {
                e0[kind][pos] -= alpha * gradient[kind][pos];
            }
        }
        for kind in 0..2 {
            for pos in 0..NUM_BITS {
                if e0[kind][pos] < 0. {
                    e0[kind][pos] = -e0[kind][pos];
                }
            }
        }
//        println!("{:?}", &e0);
        let err = (gradient.iter().flatten().map(|g: &f32| g.powi(2)).sum::<f32>() / (gradient.len() as f32)).sqrt();
//        let err = objective(x, y, &e0, max_hamming_distance, x_ind);
        println!("{}: {:?}\t{:?}", i, err, alpha);
        dbg!(&e0);
        if i >= min_iter {
            if err < target_error {
                return (*e0, err, i);
            }
        }
    }
    return (*e0, std::f32::INFINITY, max_iter);
}

pub fn estimate_errors(x: ExprV, y: ExprV, e: &Errors, max_hamming_distance: usize) -> Result<Errors, ()> {
    let x_ind: Vec<usize> = x.iter().enumerate().filter(|&(i, &v)| v != 0.).map(|(i, _)| i).collect();
    let max_iter = 256;
    let mut e0 = [[0.; NUM_BITS], [0.; NUM_BITS]];
    e0.clone_from_slice(e);
    let (e, _error, num_iters) = _gradient_descent_hardcoded(&mut e0, y, x, 1e-8, 1., max_hamming_distance, 4, max_iter, &x_ind[..]);
    if num_iters == max_iter {
        Result::Err(())
    } else {
        Result::Ok(e)
    }
}

pub fn estimate_expression(mat: &CSR, x_est: ExprV, y: ExprV, max_hamming_distance: usize, keep_zeros: bool) -> Result<Expr, ()> {
    if let Ok((x, _it, _err)) = csr_successive_overrelaxation(mat, y, x_est, 1.25, 1e-7, 256, keep_zeros) {
        Ok(x)
    } else {
        Err(())
    }
}

arg_enum! {
    #[derive(Debug, Clone)]
    pub enum Mode {
        ErrorsThenExpression,
        ExpressionThenErrors,
    }
}
