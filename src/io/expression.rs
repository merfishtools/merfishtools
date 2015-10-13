#[derive(RustcEncodable, RustcDecodable)]
pub struct Record {
    pub experiment: u32,
    pub cell: u32,
    pub feature: String
}
