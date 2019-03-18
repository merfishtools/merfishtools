use crate::model::la::common::{Errors, NUM_BITS, NUM_CODES};
use crate::model::la::hamming::hamming_distance16;
use crate::model::la::problem::prob;
use ndarray::Array1;
use rayon::prelude::*;
use std::ops::Mul;

// NNZ = {i: num_entries(i) for i in range(2, 16 + 1)}
// Pre-calculate the number of nonzero entries for CSR representation
// This is basically the result of
// `Counter(hamming_distance(i, j) for i, j in product(range(NUM_CODES), range(NUM_CODES)))`
// (with arrays instead of dicts)
const NNZ: [[usize; 17]; 17] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // 0
    [2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // 1
    [4, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // 2
    [8, 24, 24, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  // ...
    [16, 64, 96, 64, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [32, 160, 320, 320, 160, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [64, 384, 960, 1280, 960, 384, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [128, 896, 2688, 4480, 4480, 2688, 896, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [256, 2048, 7168, 14336, 17920, 14336, 7168, 2048, 256, 0, 0, 0, 0, 0, 0, 0, 0],
    [512, 4608, 18432, 43008, 64512, 64512, 43008, 18432, 4608, 512, 0, 0, 0, 0, 0, 0, 0],
    [1024, 10240, 46080, 122880, 215040, 258048, 215040, 122880, 46080, 10240, 1024, 0, 0, 0, 0, 0, 0],
    [2048, 22528, 112640, 337920, 675840, 946176, 946176, 675840, 337920, 112640, 22528, 2048, 0, 0, 0, 0, 0],
    [4096, 49152, 270336, 901120, 2027520, 3244032, 3784704, 3244032, 2027520, 901120, 270336, 49152, 4096, 0, 0, 0, 0],
    [8192, 106496, 638976, 2342912, 5857280, 10543104, 14057472, 14057472, 10543104, 5857280, 2342912, 638976, 106496, 8192, 0, 0, 0],
    [16384, 229376, 1490944, 5963776, 16400384, 32800768, 49201152, 56229888, 49201152, 32800768, 16400384, 5963776, 1490944, 229376, 16384, 0, 0],
    [32768, 491520, 3440640, 14909440, 44728320, 98402304, 164003840, 210862080, 210862080, 164003840, 98402304, 44728320, 14909440, 3440640, 491520, 32768, 0],
    [65536, 1048576, 7864320, 36700160, 119275520, 286261248, 524812288, 749731840, 843448320, 749731840, 524812288, 286261248, 119275520, 36700160, 7864320, 1048576, 65536]
];

pub struct CSR {
    data: Array1<f32>,
    columns: Array1<usize>,
    indptr: Array1<usize>,
}

impl CSR {
//    fn with_capacity(nnz: usize, n: usize) -> CSR {
//        CSR {
//            data: Array1::zeros(nnz),
//            columns: Array1::zeros(nnz),
//            indptr: Array1::zeros(n + 1),
//        }
//    }

    fn from(data: Vec<f32>, columns: Vec<usize>, indptr: Vec<usize>) -> CSR {
        CSR {
            data: Array1::from(data),
            columns: Array1::from(columns),
            indptr: Array1::from(indptr),
        }
    }
}

impl Mul<&[f32]> for CSR {
    type Output = Vec<f32>;

    fn mul(self, rhs: &[f32]) -> Self::Output {
        let mut result: Vec<f32> = vec![0.; rhs.len()];
        result
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, v)| {
                let x = (self.indptr[i]..self.indptr[i + 1]).into_par_iter()
                    .map(|k| self.data[k] * rhs[self.columns[k]])
                    .sum::<f32>();
                *v += x;
            });
        result
    }
}

impl Mul<&[f32]> for &CSR {
    type Output = Vec<f32>;

    fn mul(self, rhs: &[f32]) -> Self::Output {
        let mut result: Vec<f32> = vec![0.; rhs.len()];
        result
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, v)| {
                let x = (self.indptr[i]..self.indptr[i + 1]).into_par_iter()
                    .map(|k| self.data[k] * rhs[self.columns[k]])
                    .sum::<f32>();
                *v += x;
            });
        result
    }
}

pub fn csr_error_matrix(e: &Errors, max_hamming_distance: usize) -> CSR {
    // Build a sparse representation of the transition matrix in CSR format. Since no element is actually nonzero,
    // a threshold at which to consider an entry small enough to be negligible needs to be specified. This is achieved
    // by setting `max_hamming_distance` to an appropriate value , since matrix entries with large hamming_distance
    // correspond to small probabilities.
    let hamming_dist_count = NNZ[NUM_BITS];  // [count_entries_where(hamming_dist == i) for i in range(0, NUM_BITS + 1)]
    let nnz: usize = hamming_dist_count[..max_hamming_distance + 1].iter().sum(); // number of nonzero entries
    let mut data = vec![0.; nnz];
    let mut indices = vec![0; nnz];
    let mut indptr = vec![0; NUM_CODES + 1];
    indptr[0] = 0;
    let mut global_inc = 0;
    (0..NUM_CODES).for_each(|i| {
        let mut col_inc = 0;
        (0..NUM_CODES).filter(|&j| hamming_distance16(i, j) <= max_hamming_distance)
            .for_each(|j| {
                data[global_inc] = prob(j, i, e);
                indices[global_inc] = j;
                col_inc += 1;
                global_inc += 1;
            });
        indptr[i + 1] = indptr[i] + col_inc
    });
    CSR::from(data, indices, indptr)
}

pub fn error_dot(e: &Errors, y: &[f32], max_hamming_distance: usize) -> Vec<f32> {
    let mut data = vec![0.; NUM_CODES];
    data.par_iter_mut()
        .enumerate()
        .for_each(|(i, v)| {
            (0..NUM_CODES).filter(|j| hamming_distance16(i, *j) <= max_hamming_distance)
                .for_each(|j| *v += prob(j, i, e) * y[j]);
        });
    data
}

pub fn rmse(a: &[f32], b: &[f32]) -> f32 {
    assert_eq!(a.len(), b.len());
    let sqsum: f32 = a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x - y).powi(2)).sum();
    (sqsum / (a.len() as f32)).sqrt()
}

pub fn csr_successive_over_relaxation(a: &CSR,
                                      x: &mut [f32],
                                      y: &[f32],
                                      w: f32,
                                      eps: f32,
                                      max_iter: usize) -> Result<(Vec<f32>, usize, f32), (usize, f32)> {
    assert!(w > 0.);
    assert!(w < 2.);
    assert_eq!(x.len(), y.len());
    let mut error = rmse(&(a * x), y);
    let nrows = y.len();
    for it in 0..max_iter {
        for row in 0..nrows {
            let mut sigma = 0.;
            for k in a.indptr[row]..a.indptr[row + 1] {
                let col = a.columns[k];
                if row != col {
                    sigma += a.data[k] * x[col];
                }
            }
            let mut diag = 1.;
            for k in a.indptr[row]..a.indptr[row + 1] {
                if row == a.columns[k] {
                    diag = a.data[k]
                }
            }
            x[row] = ((1. - w) * x[row] + (w / diag) * (y[row] - sigma));
        }
        error = rmse(&(a * x), y);
        if error < eps {
            return Ok((x.to_vec(), it, error));
        }
    }
    Err((max_iter, error))
}

pub fn error_successive_overrelaxation(e: &Errors,
                                       max_hamming_distance: usize,
                                       y: &[f32],
                                       w: f32,
                                       eps: f32,
                                       max_iter: usize) -> Result<(Vec<f32>, usize, f32), (usize, f32)> {
    assert!(w > 0.);
    assert!(w < 2.);
    assert_eq!(y.len(), NUM_CODES);
    let mut x = &mut [0.; NUM_CODES];
    x.copy_from_slice(y);
    let mut error = rmse(&error_dot(e, x, max_hamming_distance), y);
    let nrows = y.len();
    for it in 0..max_iter {
        for row in 0..nrows {
            let sigma: f32 = (0..NUM_CODES)
                .filter(|&k| row != k && hamming_distance16(row, k) <= max_hamming_distance)
                .map(|col| prob(row, col, e) * x[col])
                .sum();
            let diag = prob(row, row, e);
            x[row] = ((1. - w) * x[row] + (w / diag) * (y[row] - sigma));
        }
        error = rmse(&error_dot(e, x, max_hamming_distance), y);
        if error < eps {
            return Ok((x.to_vec(), it, error));
        }
    }
    Err((max_iter, error))
}