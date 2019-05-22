use rayon::prelude::*;

use crate::model::la::common::{hamming_distance, Errors, ExprV};

fn _prob(i: usize, j: usize, k: usize, e: &Errors) -> f32 {
    let i_k = (i >> k) & 1; // kth bit in barcode i
    let j_k = (j >> k) & 1; // kth bit in barcode j
    let mut p = e[[i_k, k]]; // select relevant positional error probability
    let xor = i_k ^ j_k; // a_ijk can be written as (sign * e[i(k), k] + add) where sign in {-1, +1} and add in {0, 1}
    let add = ((!xor) & 1) as f32;
    let sign = (xor << 1) as f32 - 1.;
    p *= sign;
    p += add;
    p
}

pub fn prob(i: usize, j: usize, e: &Errors, num_bits: usize) -> f32 {
    (0..num_bits).map(|k| _prob(i, j, k, e)).product()
}

fn dprob_inf(i: usize, j: usize, e: &Errors, pos: usize, kind: usize) -> f32 {
    let i_k = (i >> pos) & 1;
    if i_k != kind {
        // we return infinity here so we can write p / dp directly
        // and thus do not need to check for dp == 0. first
        return std::f32::INFINITY;
    }
    let j_k = (j >> pos) & 1;
    let dsign = ((i_k ^ j_k) * 2) as f32 - 1.;
    let p = _prob(i, j, pos, e);
    p * dsign
}

pub fn objective(
    x: ExprV,
    y: ExprV,
    e: &Errors,
    max_hamming_distance: usize,
    x_ind: &[usize],
    y_ind: &[usize],
    num_bits: usize
) -> f32 {
    let r: f32 = y_ind
        .into_par_iter()
        .cloned()
        .map(|k| {
            let s: f32 = x_ind
                .iter()
                .cloned()
                .filter(|&j| hamming_distance(j, k) <= max_hamming_distance)
                .map(|j| prob(j, k, e, num_bits) * x[j])
                .sum();
            (y[k] - s).powi(2)
        })
        .sum();
    r
}

pub fn partial_objective(
    x: ExprV,
    y: ExprV,
    e: &Errors,
    pos: usize,
    kind: usize,
    max_hamming_distance: usize,
    x_ind: &[usize],
    y_ind: &[usize],
    num_bits: usize,
) -> f32 {
    let r: f32 = y_ind
        .into_par_iter()
        .cloned()
        .map(|k| {
            let (f, g) = x_ind
                .iter()
                .cloned()
                .filter(|&j| hamming_distance(j, k) <= max_hamming_distance)
                .map(|j| {
                    let p = prob(j, k, e, num_bits);
                    let dp = dprob_inf(j, k, e, pos, kind);
                    (p * x[j], p / dp * x[j])
                })
                // .reduce(|| (0., 0.), |(x0, y0), (x1, y1)| (x0 + x1, y0 + y1));  // par_iter version of fold below
                .fold((0., 0.), |(x0, y0), (x1, y1)| (x0 + x1, y0 + y1));
            g * (y[k] - f)
        })
        .sum();
    -r
}
