use std::cmp::Ordering;

use ndarray::prelude::{Array1, ArrayView1, ArrayViewMut1};
use ndarray::{Array2, ArrayView2, ArrayViewMut2};

pub type Errors = Array2<f32>;
// [[f32; NUM_BITS]; 2];
pub type ErrorsV<'a> = ArrayView2<'a, f32>;
pub type ErrorsVMut<'a> = ArrayViewMut2<'a, f32>;
pub type ExprV<'a> = ArrayView1<'a, f32>;
pub type ExprVMut<'a> = &'a mut ArrayViewMut1<'a, f32>;
pub type Expr = Array1<f32>;

pub fn hamming_distance(a: usize, b: usize) -> usize {
    let x = a ^ b;
    x.count_ones() as usize
}

pub fn median(x: &[f32]) -> Option<f32> {
    let mut x = x.to_owned();
    x.sort_by(|&a, &b| match a - b {
        x if x < 0. => Ordering::Less,
        x if x > 0. => Ordering::Greater,
        _ => Ordering::Equal,
    });
    Some(match x.len() {
        0 => return None,
        1 => x[0],
        len if len % 2 == 0 => {
            let x1 = x[len / 2 - 1];
            let x2 = x[len / 2];
            (x1 + x2) / 2.
        }
        len => x[len / 2],
    })
}
