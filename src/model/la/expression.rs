use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;

use clap::arg_enum;
use failure::Error;
use ndarray::prelude::*;
use rand;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rayon::prelude::*;
use regex::Regex;

use derive_new::new;
use num_traits::float::Float;

use crate::cli::Expression;
use crate::io::merfishdata;
use crate::io::merfishdata::MerfishRecord;
use crate::io::simple_codebook::SimpleCodebook;
use crate::model::la::common::{hamming_distance, Errors, Expr, ExprV};
use crate::model::la::matrix::{csr_error_matrix, csr_successive_overrelaxation, CSR};
use crate::model::la::problem::{objective, partial_objective};
use crate::simulation::binary;

#[derive(Serialize, Deserialize, new)]
pub struct RecordEstimate {
    cell: String,
    #[serde(with = "binary")]
    barcode: u16,
    #[serde(rename = "count")]
    expression: f32,
}

#[derive(new)]
pub struct ExpressionT {
    p0: Vec<f64>,
    p1: Vec<f64>,
    codebook_path: String,
    estimate_path: Option<String>,
    errors_path: Option<String>,
    cells: Regex,
    max_hamming_distance: usize,
    omega: f32,
    mode: Mode,
    seed: u64,
    num_bits: usize,
    #[new(default)]
    corrected_counts: HashMap<String, HashMap<u16, usize>>,
    #[new(default)]
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
        let codebook = SimpleCodebook::from_file(&self.codebook_path()).unwrap();
        let num_bits = self.num_bits;
        let num_codes = 1 << num_bits;
        let mut unexpressed_codewords: HashSet<u16> = codebook
            .records()
            .iter()
            .filter(|&r| !r.expressed())
            .map(|r| r.codeword())
            .collect();
        // make sure all zeros and all ones are always considered
        unexpressed_codewords.insert(0);
        unexpressed_codewords.insert((num_codes - 1) as u16);
        let y_ind: Vec<usize> = unexpressed_codewords.iter().map(|&v| v as usize).collect();

        let expressed_codewords: HashSet<usize> = codebook
            .records()
            .iter()
            .filter(|r| r.expressed())
            .map(|r| r.codeword() as usize)
            .collect();
        let _unexpressed_codewords: HashSet<usize> = codebook
            .records()
            .iter()
            .filter(|r| !r.expressed())
            .map(|r| r.codeword() as usize)
            .collect();
        let mut unused_codewords: HashSet<usize> = (0..num_codes).collect();
        for cw in &expressed_codewords {
            unused_codewords.remove(cw);
        }
        let non_codewords = Vec::from_iter(unused_codewords.iter().cloned());
        //        let non_codewords = Vec::from_iter(unexpressed_codewords.iter().cloned());
        let raw_counts = &self.raw_counts;
        let corrected_counts = &self.corrected_counts;
        let max_hamming_distance = self.max_hamming_distance;
        let mode = &self.mode;

        let mut e: Errors = Array2::from_shape_fn((2, num_bits), |(kind, pos)| {
            match kind {
                0 => self.p0[pos] as f32,
                1 => self.p1[pos] as f32,
                _ => panic!(),
            }
        });

        let mut x_est: Expr = Array1::<f32>::zeros(num_codes);
        let mut y: Expr = Array1::<f32>::zeros(num_codes);
        for (_cell, feature_counts) in corrected_counts {
            for (&barcode, &count) in feature_counts {
                x_est[barcode as usize] += count as f32;
            }
        }
        //        let uniform = rand::distributions::Uniform::new(0f32, 1.);
        //        let mut rng = StdRng::seed_from_u64(self.seed);
        for (_cell, feature_counts) in raw_counts {
            for (&barcode, &count) in feature_counts {
                y[barcode as usize] += count as f32; // + rng.sample(&uniform);
            }
        }

        let _magnitude_x = x_est.sum();
        adjust(&mut x_est, &non_codewords);
        let magnitude_y = y.sum();
        y /= magnitude_y;

        let yv = y.view();
        let mut x: Expr = x_est;

        let w = self.omega;

        match mode {
            Mode::ErrorsThenExpression => {
                for hamming_dist in 1..=max_hamming_distance {
                    println!("\n=== {:?} ===", hamming_dist);
                    let num_repeat_iters = if hamming_dist == max_hamming_distance {
                        3
                    } else {
                        1
                    };
                    for _ in 0..num_repeat_iters {
                        // set < 0. to 0., set non_codewords to 0., normalize to 1
                        adjust(&mut x, &non_codewords);

                        println!("Estimating positional error probabilities via GD");
                        e = estimate_errors(x.view(), yv, &y_ind[..], &e, hamming_dist, num_bits)
                            .expect("Failed estimating errors");
                        dbg!(&e);

                        println!("Constructing CSR transition matrix");
                        let mat = csr_error_matrix(&e, hamming_dist.max(2), num_bits);

                        println!("Estimating true expression via SOR");
                        x = estimate_expression(&mat, x.view(), yv, w, false)
                            .expect("Failed estimating expression");
                    }
                }
            }
            Mode::ExpressionThenErrors => {
                for hamming_dist in 1..=max_hamming_distance {
                    println!("\n=== {:?} ===", hamming_dist);
                    let num_repeat_iters = if hamming_dist == max_hamming_distance {
                        3
                    } else {
                        1
                    };
                    for _ in 0..num_repeat_iters {
                        println!("Constructing CSR transition matrix");
                        let mat = csr_error_matrix(&e, hamming_dist.max(2), num_bits);

                        println!("Estimating true expression via SOR");
                        adjust(&mut x, &non_codewords);
                        x = estimate_expression(&mat, x.view(), yv, w, false)
                            .expect("Failed estimating expression");
                        adjust(&mut x, &non_codewords);

                        println!("Estimating positional error probabilities via GD");
                        e = estimate_errors(x.view(), yv, &y_ind[..], &e, hamming_dist, num_bits)
                            .expect("Failed estimating errors");
                        dbg!(&e);
                    }
                }
            }
        };
        let mat = csr_error_matrix(&e, 4, num_bits);
        let estimate_writer = self.estimate_path.as_ref().map(|path| {
            csv::WriterBuilder::new()
                .delimiter(b'\t')
                .from_path(path)
                .unwrap()
        });
        let all_records: Vec<RecordEstimate> = self
            .raw_counts
            .par_iter()
            .flat_map(|(cell, raw_countmap)| {
                let corrected_countmap = self.corrected_counts.get(cell).unwrap();
                let mut y: Expr = Expr::zeros((num_codes, ));
                raw_countmap.iter().for_each(|(&barcode, &count)| {
                    y[barcode as usize] += count as f32;
                });
                let mut x: Expr = Expr::zeros((num_codes, ));
                corrected_countmap.iter().for_each(|(&barcode, &count)| {
                    x[barcode as usize] += count as f32;
                });

                let magnitude_x = x.sum();
                adjust(&mut x, &non_codewords);
                let magnitude_y = y.sum();
                //                dbg!(magnitude_x);
                //                dbg!(magnitude_y);
                y /= magnitude_y;
                x = estimate_expression(&mat, x.view(), y.view(), w, false)
                    .expect("Failed estimating expression");
                adjust(&mut x, &non_codewords);
                x *= magnitude_x;
                x.iter()
                    .enumerate()
                    .filter(|(_, &v)| v > 0.)
                    .map(|(i, &v)| RecordEstimate::new(cell.to_owned(), i as u16, v))
                    .collect::<Vec<RecordEstimate>>()
            })
            .collect();
        if let Some(mut writer) = estimate_writer {
            for r in all_records {
                writer.serialize(r)?;
            }
        }
        if let Some(errors_path) = &self.errors_path {
            std::fs::write(errors_path, format!("{:?}", &e))?;
        }

        Ok(())
    }
}

fn adjust(x: &mut Expr, non_codewords: &[usize]) {
    x.iter_mut().filter(|&&mut v| v < 0.).for_each(|v| *v = 0.);
    for &cw in non_codewords {
        x[cw] = 0.;
    }
    *x /= x.sum();
}

impl ExpressionT {
    pub fn load_counts<'a, R>(
        &mut self,
        reader: &'a mut R,
        format: merfishdata::Format,
    ) -> Result<(), Error>
        where
            R: crate::io::merfishdata::Reader<'a>,
    {
        let codebook =
            &crate::io::simple_codebook::SimpleCodebook::from_file(&self.codebook_path())?;

        let mut corrected_counts: HashMap<String, HashMap<u16, usize>> = HashMap::new();
        let mut raw_counts: HashMap<String, HashMap<u16, usize>> = HashMap::new();

        let mut rng = StdRng::seed_from_u64(self.seed);

        let expressed_codewords: Vec<u16> = codebook
            .records()
            .iter()
            .filter(|&r| r.expressed())
            .map(|r| r.codeword())
            .collect();
        let _unexpressed_codewords: Vec<u16> = codebook
            .records()
            .iter()
            .filter(|&r| !r.expressed())
            .map(|r| r.codeword())
            .collect();

        for rec in reader.records() {
            let record = rec?;
            //            if !self.cells().is_match(&record.cell_name()) && codebook.contains(&record.feature_name()) {
            //                continue;
            //            }
            //            if !codebook.record(record.feature_id() as usize).expressed() {
            //                continue;
            //            }
            let barcode = record.barcode(None);
            let mut raw_barcode = barcode;
            let mut uncorrected_barcode = barcode;
            match format {
                // simulated data is perturbed ground truth data,
                // i.e. no error correction has been performed
                // ... so we do that now, if at all possible
                merfishdata::Format::Simulation => {
                    for dist in 1..=4u8 {
                        if record.hamming_dist() == dist {
                            let closest: Vec<u16> = expressed_codewords
                                .iter()
                                .map(|&cw| {
                                    (cw, hamming_distance(cw as usize, raw_barcode as usize))
                                })
                                .filter(|(_cw, d)| *d == dist as usize)
                                .map(|(cw, _)| cw)
                                .collect();
                            if closest.len() == 1 {
                                dbg!("FOUND ONE!");
                                raw_barcode = closest[0];
                                break;
                            }
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
                                let weights: Vec<_> = (0..self.num_bits)
                                    .map(|i| {
                                        let bit = (barcode >> i) & 1;
                                        let p = match bit {
                                            0 => self.p0[i],
                                            1 => self.p1[i],
                                            _ => panic!("bug: bit can only be 0 or 1"),
                                        };
                                        p
                                    })
                                    .collect();
                                let wc = WeightedIndex::new(weights).unwrap();
                                rng.sample(&wc) as u8
                            }
                        };
                        uncorrected_barcode ^= 1 << error_bit;
                    }
                }
            }

            let count = record.count();
            let cell_counts_corrected = corrected_counts
                .entry(record.cell_name())
                .or_insert_with(HashMap::new);
            let feature_counts_corrected = cell_counts_corrected.entry(raw_barcode).or_insert(0);
            *feature_counts_corrected += count;

            let cell_counts_raw = raw_counts
                .entry(record.cell_name())
                .or_insert_with(HashMap::new);
            let feature_counts_raw = cell_counts_raw.entry(uncorrected_barcode).or_insert(0);
            *feature_counts_raw += count;
        }
        self.corrected_counts = corrected_counts;
        self.raw_counts = raw_counts;
        Ok(())
    }
}

fn _gradient_descent_hardcoded(
    e: Errors,
    y: ExprV,
    x: ExprV,
    target_error: f32,
    alpha: f32,
    max_hamming_distance: usize,
    min_iter: usize,
    max_iter: usize,
    x_ind: &[usize],
    y_ind: &[usize],
    num_bits: usize,
) -> (Errors, f32, usize) {
    let mut gradient: Errors;
    let mut e0 = e;
    let mult: f32 = 0.5;
    let c1: f32 = 1e-4;
    let mut last_err = f32::infinity();
    let alpha1 = alpha;
    for i in 0..max_iter {
        //        println!("Calculating gradient...");
        gradient = Array2::from_shape_fn((2, num_bits), |(kind, pos)| {
            partial_objective(x, y, &e0, pos, kind, max_hamming_distance, x_ind, y_ind, num_bits)
        });
        //        dbg!(&gradient);
        //        println!("Estimating stepsize...");
        let o2 = objective(x, y, &e0, max_hamming_distance, x_ind, y_ind, num_bits);
        // gradv → descent direction: -partial_objective^T · partial_objective
        let gradv = -gradient.map(|v| v.powi(2)).sum();

        //        alpha1 = 0.5 * alpha1 + 0.5 * (0.001 / (gradient.iter().flatten().map(|v: &f32| v.abs()).sum::<f32>() / (2 * NUM_BITS) as f32));
        //        dbg!(&gradient);
        //        dbg!((o2, gradv));
        let mut alpha = alpha1;
        for _ in 0..17 {
            let mut e1 = e0.to_owned();
            e1 = e1 - alpha * &gradient;
            e1.mapv_inplace(f32::abs);
            // e1 = e0 + alpha · gradv
            let o1 = objective(x, y, &e1, max_hamming_distance, x_ind, y_ind, num_bits);
            let o2 = o2 + c1 * alpha * gradv;
            // armijo rule →  f(e0 + alpha · gradv) <= f(e0) + c1 · alpha · gradv?
            if o1 <= o2 {
                break;
            }
            alpha *= mult;
        }
        e0 = e0 - alpha * &gradient;

        //        let magnitudes: Vec<f32> = e0[0].iter().map(|v| v.log10().ceil()).collect();
        //        let large_magnitudes: Vec<f32> = e0[0].iter().zip(magnitudes.iter()).filter(|(_, &m)| m > -2.).map(|(&p, _)| p).collect();
        //        let median_e0 = median(&large_magnitudes[..]).unwrap();
        //        for p0 in e0[0].iter_mut() {
        //            let m = p0.log10().ceil();
        //            if m <= -2. {
        //                *p0 = median_e0;
        //            }
        //        }
        //
        //        let magnitudes: Vec<f32> = e0[1].iter().map(|v| v.log10().ceil()).collect();
        //        let large_magnitudes: Vec<f32> = e0[1].iter().zip(magnitudes.iter()).filter(|(_, &m)| m > -2.).map(|(&p, _)| p).collect();
        //        let median_e0 = median(&large_magnitudes[..]).unwrap();
        //        for p0 in e0[1].iter_mut() {
        //            let m = p0.log10().ceil();
        //            if m <= -2. {
        //                *p0 = median_e0;
        //            }
        //        }

        e0.mapv_inplace(|v| v.abs());
        e0.mapv_inplace(|v| if v > 0.5 { 0.5 } else { v });

        //        println!("{:?}", &e0);
        let err = (gradient.map(|v| v.powi(2)).sum() / (gradient.len() as f32)).sqrt();
        //        let err = objective(x, y, &e0, max_hamming_distance, x_ind);
        //        println!("{}: {:e}\t{:?}", i, err, alpha);
        //        dbg!(&e0);
        if i >= min_iter {
            if last_err < err {
                return (e0.to_owned(), err, i);
            } else {
                last_err = err;
            }
            if err < target_error {
                return (e0.to_owned(), err, i);
            }
        }
    }
    return (e0.to_owned(), std::f32::INFINITY, max_iter);
}

pub fn estimate_errors(
    x: ExprV,
    y: ExprV,
    y_ind: &[usize],
    e: &Errors,
    max_hamming_distance: usize,
    num_bits: usize,
) -> Result<Errors, ()> {
    let x_ind: Vec<usize> = x
        .iter()
        .enumerate()
        .filter(|&(_i, &v)| v != 0.)
        .map(|(i, _)| i)
        .collect();
    let max_iter = 256;
    let e0 = e.to_owned();
    let (e, _error, num_iters) = _gradient_descent_hardcoded(
        e0,
        y,
        x,
        1e-9,
        5e5,
        max_hamming_distance,
        4,
        max_iter,
        &x_ind[..],
        &y_ind[..],
        num_bits,
    );
    if num_iters == max_iter {
        Result::Err(())
    } else {
        Result::Ok(e)
    }
}

pub fn estimate_expression(
    mat: &CSR,
    x_est: ExprV,
    y: ExprV,
    w: f32,
    keep_zeros: bool,
) -> Result<Expr, ()> {
    if let Ok((x, _it, _err)) =
    csr_successive_overrelaxation(mat, y, x_est, w, 1e-7, 256, keep_zeros)
    {
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
