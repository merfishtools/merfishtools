use std::io;
use std::fs;
use std::path::Path;
use std::collections::HashMap;
use std::hash::Hash;

use itertools::Itertools;
use csv;
use rustc_serialize::{Decodable, Encodable};

use bio::stats::logprobs::LogProb;


#[derive(RustcEncodable, RustcDecodable)]
pub struct Record<F: Decodable + Encodable, V: Decodable + Encodable> {
    pub feature: F,
    pub value: V,
    pub prob: LogProb
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

    pub fn records<'a, F: Decodable + Encodable, V: Decodable + Encodable>(&'a mut self) -> csv::DecodedRecords<'a, R, Record<F, V>> {
        self.inner.decode()
    }

    pub fn pmf<F: Decodable + Encodable + Eq + Hash + Clone, V: Decodable + Encodable + Clone>(&mut self) -> HashMap<F, Vec<(V, LogProb)>> {
        let mut pmf = HashMap::new();
        for (feature, records) in self.records::<F, V>().map(|res| res.ok().expect("Error reading record")).group_by(
            |rec| rec.feature.clone()
        ) {
            let items = records.iter().map(|rec| (rec.value.clone(), rec.prob.clone())).collect_vec();
            pmf.insert(feature, items);
        }
        pmf
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

    pub fn write_header(&mut self, header: &[&str]) -> csv::Result<()> {
        self.inner.write(header.into_iter())
    }

    pub fn write<F: Decodable + Encodable, V: Decodable + Encodable>(&mut self, record: &Record<F, V>) -> csv::Result<()> {
        self.inner.encode(record)
    }
}
