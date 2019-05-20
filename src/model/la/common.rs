use ndarray::prelude::{Array1, ArrayView1, ArrayViewMut1};

#[cfg(feature = "14bit")]
pub const NUM_BITS: usize = 14;

#[cfg(not(feature = "14bit"))]
pub const NUM_BITS: usize = 16;
pub const NUM_CODES: usize = 2 << (NUM_BITS - 1);

pub type Errors = [[f32; NUM_BITS]; 2];
pub type ExprV<'a> = ArrayView1<'a, f32>;
pub type ExprVMut<'a> = &'a mut ArrayViewMut1<'a, f32>;
pub type Expr = Array1<f32>;

pub fn hamming_distance(a: usize, b: usize) -> usize {
    let x = a ^ b;
    x.count_ones() as usize
}
