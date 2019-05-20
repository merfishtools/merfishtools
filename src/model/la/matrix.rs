use std::ops::Mul;

use ndarray::{Array1, ArrayView1, Axis};
use ndarray::prelude::*;
use ndarray_parallel::*;
use rayon::prelude::*;

use crate::model::la::common::{Errors, Expr, ExprV, NUM_BITS, NUM_CODES, hamming_distance};
use crate::model::la::problem::prob;

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

macro_rules! impl_mul_vec {
   ($base_type: ty, $in_type: ty, $out_type_inner: ty) => {
        impl<'a> Mul<$in_type> for $base_type {
            type Output = Vec<$out_type_inner>;

            fn mul(self, rhs: $in_type) -> Self::Output {
                let mut result: Vec<$out_type_inner> = vec![0.; rhs.len()];
                result
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(i, v)| {
                        let x = (self.indptr[i]..self.indptr[i + 1]).into_par_iter()
                            .map(|k| self.data[k] * rhs[self.columns[k]])
                            .sum::<$out_type_inner>();
                        *v += x;
                    });
                result
            }
        }
    }
}

macro_rules! impl_mul_arr {
    ($base_type: ty, $in_type: ty, $out_type_inner: ty) => {
        impl<'a> Mul<$in_type> for $base_type {
            type Output = Array1<$out_type_inner>;

            fn mul(self, rhs: $in_type) -> Self::Output {
                let mut result: Array1<$out_type_inner> = Array1::<$out_type_inner>::zeros((rhs.len(),));
                result
                    .axis_iter_mut(Axis(0))
                    .into_par_iter()
                    .enumerate()
                    .for_each(|(i, mut v)| {
                        let x = (self.indptr[i]..self.indptr[i + 1]).into_par_iter()
                            .map(|k| self.data[k] * rhs[self.columns[k]])
                            .sum::<f32>();
                        // v is actually a single value
                        v.mapv_inplace(|_s| x);
                });
                result
            }
        }
    }
}

impl_mul_vec!(CSR, &[f32], f32);
impl_mul_vec!(&CSR, &[f32], f32);
impl_mul_arr!(CSR, ArrayView1<'a, f32>, f32);
impl_mul_arr!(&CSR, ArrayView1<'a, f32>, f32);
impl_mul_arr!(CSR, ArrayViewMut1<'a, f32>, f32);
impl_mul_arr!(&CSR, ArrayViewMut1<'a, f32>, f32);
impl_mul_arr!(CSR, &ArrayView1<'a, f32>, f32);
impl_mul_arr!(&CSR, &ArrayView1<'a, f32>, f32);
impl_mul_arr!(CSR, &mut ArrayViewMut1<'a, f32>, f32);
impl_mul_arr!(&CSR, &mut ArrayViewMut1<'a, f32>, f32);

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
        (0..NUM_CODES).filter(|&j| hamming_distance(i, j) <= max_hamming_distance)
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

pub fn error_dot(e: &Errors, y: ArrayView1<f32>, max_hamming_distance: usize) -> Array1<f32> {
    let mut data = Array1::zeros((y.len(), ));
    data.axis_iter_mut(Axis(0))
        .into_par_iter()
        .enumerate()
        .for_each(|(i, mut v)| {
            (0..NUM_CODES).filter(|j| hamming_distance(i, *j) <= max_hamming_distance)
                .for_each(|j| {
                    v.mapv_inplace(|_| prob(j, i, e) * y[j]);
                });
        });
    data
}

pub fn rmse(a: ArrayView1<f32>, b: ArrayView1<f32>) -> f32 {
    assert_eq!(a.len(), b.len());
    let sqsum = (&a - &b).map(|v| v.powi(2)).scalar_sum();
    (sqsum / (a.len() as f32)).sqrt()
}

pub fn csr_successive_overrelaxation(a: &CSR,
                                     y: ExprV,
                                     x_est: ExprV,
                                     w: f32,
                                     eps: f32,
                                     max_iter: usize,
                                     keep_zeros: bool) -> Result<(Expr, usize, f32), (usize, f32)> {
    assert!(w > 0.);
    assert!(w < 2.);

    let mut x = x_est.to_owned();
//    let mut rng = rand::StdRng::from_seed(&[42]);
//    let std_dev = y.std_axis(Axis(0), 1.);
//    let b = std_dev.iter().cloned().take(1).collect::<Vec<f32>>()[0];
//    let mut normal = Normal::new(0., (b as f64) / 16.);
//    x.mapv_inplace(|v| (v + normal.sample(&mut rng) as f32).abs());

//    dbg!(x.slice(s![0..10; 1]));
//    dbg!(y.slice(s![0..10; 1]));

    let mut error = rmse((a * x.view()).view(), y);
    let mut last_error = std::f32::INFINITY;
//    dbg!((a * x.view()).view().slice(s![0..10; 1]));
//    dbg!(y.slice(s![0..10; 1]));
    let nonzero: Vec<usize> = if keep_zeros {
        y.iter().enumerate().filter(|(_, &v)| v != 0.).map(|(i, _)| i).collect()
    } else {
        (0..y.len()).collect()
    };
//    let num_rows = NUM_CODES;
    for it in 0..max_iter {
        for &row in &nonzero {
            let mut sigma = 0.;
            for k in a.indptr[row]..a.indptr[row + 1] {
                let col = a.columns[k];
                if row != col {
                    sigma += a.data[k] * &x[col];
                }
            }
            let mut diag = 1.;
            for k in a.indptr[row]..a.indptr[row + 1] {
                if row == a.columns[k] {
                    diag = a.data[k];
                    break;
                }
            }
            let value = (1. - w) * &x[row] + (w / diag) * (y[row] - sigma);
            if value >= 0. {
                x[row] = value;
            }
        }
        error = rmse((a * x.view()).view(), y);
//        dbg!(&error);
//        println!("{:?}, {:?}", it, &error);
        if error.is_infinite() || error.is_nan() {
            return Err((it, error));
        }
        if error < eps || (last_error - error).abs() <= 1e-7 {
//            println!("{:?}, {:?}", x.slice(s![0..10; 1]), y.slice(s![0..10; 1]));
            return Ok((x.to_owned(), it, error));
        }
        last_error = error;
    }
    Err((max_iter, error))
}

pub fn error_successive_overrelaxation(e: &Errors,
                                       max_hamming_distance: usize,
                                       y: ExprV,
                                       w: f32,
                                       eps: f32,
                                       max_iter: usize) -> Result<(Array1<f32>, usize, f32), (usize, f32)> {
    assert!(w > 0.);
    assert!(w < 2.);
    assert_eq!(y.len(), NUM_CODES);
    let mut x = y.to_owned();
    let mut error = rmse(error_dot(e, x.view(), max_hamming_distance).view(), y.view());
    let nrows = y.len();
    for it in 0..max_iter {
        for row in 0..nrows {
            let sigma: f32 = (0..NUM_CODES)
                .filter(|&k| row != k && hamming_distance(row, k) <= max_hamming_distance)
                .map(|col| prob(row, col, e) * x[col])
                .sum();
            let diag = prob(row, row, e);
            x[row] = (1. - w) * x[row] + (w / diag) * (y[row] - sigma);
        }
        error = rmse(error_dot(e, x.view(), max_hamming_distance).view(), y.view());
        if error < eps {
            return Ok((x, it, error));
        }
    }
    Err((max_iter, error))
}