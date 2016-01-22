pub mod merfishdata;
pub mod codebook;
pub mod pmf;
pub mod estimation;


#[derive(RustcEncodable, RustcDecodable, PartialEq, Eq, Clone, Copy, Hash)]
pub struct Cell {
    pub experiment: u32,
    pub cell: u32
}
