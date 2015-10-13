use std::io;
use std::fs;
use std::path::Path;

use csv;

use bio::stats::logprobs::Prob;

use io::PMF;

#[derive(RustcEncodable, RustcDecodable)]
pub struct Record {
    pub experiment: u32,
    pub cell: u32,
    pub feature: String,
    pub pmf: PMF<u32>
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


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(|f| Writer::from_writer(f))
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        }
    }

    pub fn write_header(&mut self) -> csv::Result<()> {
        self.inner.write(["Expmt", "Cell", "Feat", "Expr", "Prob"].into_iter())
    }

    pub fn write(&mut self, record: Record) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
