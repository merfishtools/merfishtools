use std::io;
use std::fs;
use std::path::Path;

use csv;


#[derive(RustcDecodable)]
pub struct Position {
    pub x: f32,
    pub y: f32
}


#[derive(RustcDecodable)]
pub struct Record {
    pub experiment: u32,
    pub codebook: u32,
    pub cell_id: u32,
    pub rna: String,
    pub exact_match: u8,
    pub corrected_match: u8,
    pub cell_position: Position,
    pub rna_position: Position
}


pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>
}


impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::from_reader(f))
    }
}


impl<R: io::Read> Reader<R> {
    pub fn from_reader(rdr: R) -> Self {
        Reader {
            inner: csv::Reader::from_reader(rdr).delimiter(b'\t')
        }
    }

    pub fn records<'a>(&'a mut self) -> csv::DecodedRecords<'a, R, Record> {
        self.inner.decode()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io;

    #[test]
    fn test_records() {
        let data = b"Experiment	Codebook	Cell_ID	Gene_Name	Exact_Match	Corrected_Match	Cell_Position_X	Cell_Position_Y	RNA_Position_X	RNA_Position_Y
1	1	0	SCUBE3	0	1	475.5	630.6	13146.86026973	25793.5656964
1	1	0	SCUBE3	0	1	475.5	630.6	13576.7356895	38396.4273422
";
        let mut reader = Reader::from_reader(io::Cursor::new(&data[..]));
        for r in reader.records() {
            match r {
                Ok(rec) => {
                    assert_eq!(rec.rna, "SCUBE3");
                    assert_eq!(rec.cell_position.x, 475.5);
                },
                Err(e) => panic!("{:?}", e)
            }
        }
    }
}
