use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;
use std::time::Instant;

use clap::arg_enum;
use derive_new::new;
use failure::Error;
use ndarray::prelude::*;
use num_traits::float::Float;
use rand;
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rayon::prelude::*;
use regex::Regex;

use crate::cli::Expression;
use crate::io::common::Barcode;
use crate::io::codebook::Codebook;
use crate::io::counts::{Counts, Records};
use crate::io::merfishdata;
use crate::io::merfishdata::MerfishRecord;
use crate::io::simple_codebook::SimpleCodebook;
use crate::model::la::common::{hamming_distance, Errors, Expr, ExprV};
use crate::model::la::matrix::{csr_error_matrix, csr_successive_overrelaxation, CSR};
use crate::model::la::problem::{objective, partial_objective};
use crate::simulation::binary;
use itertools::{Either, Itertools};
use std::path::Path;

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
    counts: HashMap<String, Counts>,
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
        let max_hamming_distance = self.max_hamming_distance;
        let mode = &self.mode;

        let mut e: Errors = Array2::from_shape_fn((2, num_bits), |(kind, pos)| match kind {
            0 => self.p0[pos] as f32,
            1 => self.p1[pos] as f32,
            _ => panic!(),
        });

        let mut x_est: Expr = Array1::<f32>::zeros(num_codes);
        let mut y: Expr = Array1::<f32>::zeros(num_codes);
        for (_cell, counts) in &self.counts {
            for (&(_codeword, readout), &count) in counts.corrected.iter() {
                x_est[*readout as usize] = count as f32;
            }
            for (&(_codeword, readout), &count) in counts.observed.iter() {
                y[*readout as usize] = count as f32;
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
                        info!("Adjusting x (x[x < 0] := 0, x[Blanks] := 0, ||x|| := 1)\n");
                        adjust(&mut x, &non_codewords);

                        info!("Estimating positional error probabilities via GD");
                        let mut now = Instant::now();
                        e = estimate_errors(x.view(), yv, &y_ind[..], &e, hamming_dist, num_bits)
                            .expect("Failed estimating errors");
                        dbg!(&e);
                        info!("Finished estimating errors in {:#?}.\n", now.elapsed());

                        info!("Constructing CSR transition matrix");
                        now = Instant::now();
                        let mat = csr_error_matrix(&e, hamming_dist.max(2), num_bits);
                        info!(
                            "Finished constructing CSR matrix in {:#?}.\n",
                            now.elapsed()
                        );

                        println!("Estimating true expression via SOR");
                        now = Instant::now();
                        x = estimate_expression(&mat, x.view(), yv, w, false)
                            .expect("Failed estimating expression");
                        info!("Finished estimating expression in {:#?}.\n", now.elapsed());
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
            .counts
            .par_iter()
            .flat_map(|(cell, c): (&String, &Counts)| {
                let mut x: Expr = Expr::zeros((num_codes, ));
                let mut y: Expr = Expr::zeros((num_codes, ));
                for (&(_codeword, readout), &count) in c.corrected.iter() {
                    x[*readout as usize] = count as f32;
                }
                for (&(_codeword, readout), &count) in c.observed.iter() {
                    y[*readout as usize] = count as f32;
                }

                let magnitude_x = x.sum();
                adjust(&mut x, &non_codewords);
                let magnitude_y = y.sum();
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

/// Adjusts expression `x` such that:
///   1. `x[(x < 0) | non_codewords] == 0`
///   2. `||x|| == 1`
///
/// # Examples
/// ```
/// # extern crate ndarray;
/// # use ndarray::{array, Array1};
/// # fn main() {
/// let mut x: Array1<f32> = array![0., -1., 2., 4., 2.];
/// let b = [3usize];
/// adjust(&mut x, &b);
/// assert_eq!(x, array![0., 0., 0.5, 0., 0.5]);
/// # }
/// # fn adjust(x: &mut Array1<f32>, non_codewords: &[usize]) {
/// #   x.iter_mut().filter(|&&mut v| v < 0.).for_each(|v| *v = 0.);
/// #   for &cw in non_codewords {
/// #       x[cw] = 0.;
/// #   }
/// #   *x /= x.sum();
///# }
/// ```
///
fn adjust(x: &mut Expr, non_codewords: &[usize]) {
    //    let sorted_x: Vec<(f32, usize)> = x.iter().cloned().zip(0..x.len())
    //        .sorted_by(|(v1, _), (v2, _)| (-v1).partial_cmp(&-v2).unwrap()).collect();
    //    let mut sum = 0.;
    //    for (v, i) in sorted_x {
    //        if sum >= 1. {
    //            dbg!((v, i));
    //            x[i] = 0.;
    //        } else {
    //            sum += v;
    //        }
    //    }

    x.iter_mut().filter(|&&mut v| v < 0.).for_each(|v| *v = 0.);
    for &cw in non_codewords {
        x[cw] = 0.;
    }
    *x /= x.sum();
}

impl ExpressionT {
    pub fn load_counts<'a, P: AsRef<Path>>(&mut self, path: P) -> Result<(), Error> {
        info!("Reading codebook");
        let codebook =
            &crate::io::simple_codebook::SimpleCodebook::from_file(&self.codebook_path())?;
        info!("Read codebook.");
        let (expressed_codewords, _unexpressed_codewords): (Vec<_>, Vec<_>) =
            codebook.records().iter().partition_map(|r| {
                if r.expressed() {
                    Either::Left(Barcode(r.codeword()))
                } else {
                    Either::Right(Barcode(r.codeword()))
                }
            });
        info!("Reading records...");
        let records = Records::from_path(path);
        info!("Finished reading records.");
        let counts = Counts::from_records(records.into_iter(), &expressed_codewords);
        info!("Finished counting.");
        self.counts = counts;
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
            partial_objective(
                x,
                y,
                &e0,
                pos,
                kind,
                max_hamming_distance,
                x_ind,
                y_ind,
                num_bits,
            )
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

/// Estimate transitional error probabilities given observation y and (estimated) true expression x.
///
/// Note that, during computation, all indices where `x == 0` may be skipped and only
/// those indices where `y` is *supposed* to be non-zero (→ `y_ind`) need to be considered.
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
        2,
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
