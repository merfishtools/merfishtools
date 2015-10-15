use std::io;
use std::fs;
use std::path::Path;

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


pub struct AggregatedRecord<F, V> {
    pub feature: F,
    pub pmf: Vec<(V, LogProb)>
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

    /*pub fn aggregated_records<'a, F: Decodable + Encodable + Eq + Clone, V: Decodable + Encodable + Clone>(&'a mut self) -> Box<Iterator<Item=AggregatedRecord<F, V>>> {
        Box::new(self.inner.decode().map(|res| res.ok().expect("Error reading record.")).group_by(|rec: Record<F, V>| rec.feature.clone()).map(|(feature, records): (F, Vec<Record<F, V>>)| {
            AggregatedRecord {
                feature: feature,
                pmf: records.iter().map(|rec: &Record<F, V>| (rec.value.clone(), rec.prob.clone())).collect_vec()
            }
        }))
    }*/
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
