pub const NUM_BITS: usize = 12;
pub const NUM_CODES: usize = 2 << (NUM_BITS - 1);

pub type Errors = [[f32; NUM_BITS]; 2];
