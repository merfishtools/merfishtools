use std::collections::HashMap;

use bio::stats::Prob;
use failure::Error;
use rand;
use rand::Rng;
use rand::StdRng;
use rayon::prelude::*;
use regex::Regex;

use crate::cli::Expression;
use crate::model::bayes::readout::Counts;
use crate::model::la::common::{Errors, NUM_BITS, NUM_CODES};
use crate::model::la::expression::Mode::ExpressionThenErrors;
use crate::model::la::matrix::{csr_error_matrix, error_dot, error_successive_overrelaxation};
use crate::model::la::problem::{objective, partial_objective};
use crate::io::codebook::Codebook;

#[derive(Builder)]
#[builder(pattern = "owned")]
pub struct ExpressionT {
    p0: Vec<Prob>,
    p1: Vec<Prob>,
    codebook_path: String,
    estimate_path: Option<String>,
    stats_path: Option<String>,
    threads: usize,
    cells: Regex,
    window_width: u32,
    seed: usize,
    #[builder(setter(skip))]
    counts: HashMap<String, HashMap<String, Counts>>,
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
        let max_hamming_distance = 4;  // TODO read from cli, add to struct.
        let mode = ExpressionThenErrors;  // TODO read from cli, add to struct.
        let mut e = [[0.; NUM_BITS], [0.; NUM_BITS]];
        e[0].iter_mut().zip(self.p0.iter().map(f32::from))
            .for_each(|(v, p)| *v = p.0);
        e[1].iter_mut().zip(self.p1.iter().map(f32::from))
            .for_each(|(v, p)| *v = p.0);
//        e.copy_from_slice(e0);
        let mut x = vec![0.; NUM_CODES];
        for hamming_dist in 1..=max_hamming_distance {
            match mode {
                Mode::ErrorsThenExpression => {
                    e = estimate_errors(&x, &y, &e).expect("Failed estimating errors");
                    x = estimate_expression(&e, &y, max_hamming_distance).expect("Failed estimating expression");
                }
                Mode::ExpressionThenErrors => {
                    x = estimate_expression(&e, &y, max_hamming_distance).expect("Failed estimating expression");
                    e = estimate_errors(&x, &y, &e).expect("Failed estimating errors");
                }
            }
        };
        unimplemented!()
    }
}

fn _gradient_descent_hardcoded(e: &mut Errors,
                               y: &[f32],
                               x: &[f32],
                               target_error: f32,
                               alpha: f32,
                               max_hamming_distance: usize,
                               min_iter: usize,
                               max_iter: usize,
                               x_ind: &[usize]) -> (Errors, f32, usize) {
    let mut last_err = std::f32::INFINITY;
    let mut gradient: Errors = [[0.; NUM_BITS]; 2];
    let e0 = e;
    let alpha0 = alpha;
    let mult: f32 = 0.5;
    let c1: f32 = 1e-2;
    for i in 0..max_iter {
        println!("Calculating gradient...");
        (0..2).for_each(|kind| {
            (0..NUM_BITS).for_each(|pos| {
                let grad_value = partial_objective(x, y, &e0, pos, kind, max_hamming_distance, x_ind);
                gradient[kind][pos] = grad_value;
            });
        });
        println!("Estimating stepsize...");
        let mut alpha = alpha0;
        for _ in 0..17 {
            let mut e1 = *e0;
            for kind in 0..2 {
                for pos in 0..NUM_BITS {
                    e1[kind][pos] -= alpha * gradient[kind][pos];
                }
            }
            let o1 = objective(x, y, &e1, max_hamming_distance, x_ind);
            let o2 = objective(x, y, &e0, max_hamming_distance, x_ind)
                - c1 * alpha * gradient.iter().flatten().map(|g: &f32| g.powi(2)).sum::<f32>();
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
        let err = (gradient.iter().flatten().map(|g: &f32| g.powi(2)).sum::<f32>() / (gradient.len() as f32)).sqrt();
        last_err = err;
        if i >= min_iter {
            if err < target_error {
                return (*e0, err, i);
            }
        }
    }
    return (*e0, std::f32::INFINITY, max_iter);
}

fn estimate_errors(x: &[f32], y: &[f32], e: &Errors) -> Result<Errors, ()> {
    let x_ind: Vec<usize> = x.iter().enumerate().filter(|&(i, &v)| v != 0.).map(|(i, _)| i).collect();
    let max_iter = 1024;
    let mut e0 = [[0.; NUM_BITS], [0.; NUM_BITS]];
    e0.clone_from_slice(e);
    let (e, error, num_iters) = _gradient_descent_hardcoded(&mut e0, y, x, 0.9, 0.25, 16, 4, max_iter, &x_ind[..]);
    if num_iters == max_iter {
        Result::Err(())
    } else {
        Result::Ok(e)
    }
}

fn estimate_expression(e: &Errors, y: &[f32], max_hamming_distance: usize) -> Result<Vec<f32>, ()> {
    if let Ok((x, it, err)) = error_successive_overrelaxation(e, max_hamming_distance, y, 0.75, 1e-7, 4096) {
        Ok(x)
    } else {
        Err(())
    }
}

pub enum Mode {
    ErrorsThenExpression,
    ExpressionThenErrors,
}

pub fn main(e0: &Errors, y: &[f32], max_hamming_distance: usize) -> Vec<f32> {
    let mode = ExpressionThenErrors;
    let mut e = [[0.; NUM_BITS], [0.; NUM_BITS]];
    e.copy_from_slice(e0);
    let mut x = vec![0.; NUM_CODES];
    for hamming_dist in 1..=max_hamming_distance {
        match mode {
            Mode::ErrorsThenExpression => {
                e = estimate_errors(&x, &y, &e).expect("Failed estimating errors");
                x = estimate_expression(&e, &y, max_hamming_distance).expect("Failed estimating expression");
            }
            Mode::ExpressionThenErrors => {
                x = estimate_expression(&e, &y, max_hamming_distance).expect("Failed estimating expression");
                e = estimate_errors(&x, &y, &e).expect("Failed estimating errors");
            }
        }
    }
    x
}
