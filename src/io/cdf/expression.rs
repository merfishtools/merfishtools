use std::io;
use std::fs;
use std::path::Path;
use std::collections;

use itertools::Itertools;
use csv;

use bio::stats::logprobs::LogProb;

use model::expression::CDF;


const HEADER: [&'static str; 4] = ["cell", "feat", "expr", "prob"];


pub struct CDFs {
    inner: collections::HashMap<String, Vec<CDF>>
}


impl CDFs {
    pub fn contains_feature(&self, feature: &str) -> bool {
        self.inner.contains_key(feature)
    }

    pub fn get(&self, feature: &str) -> Option<&Vec<CDF>> {
        self.inner.get(feature)
    }

    pub fn get_mut(&mut self, feature: &str) -> Option<&mut Vec<CDF>> {
        self.inner.get_mut(feature)
    }

    pub fn features(&self) -> collections::hash_map::Keys<String, Vec<CDF>>  {
        self.inner.keys()
    }

    pub fn iter(&self) -> collections::hash_map::Iter<String, Vec<CDF>> {
        self.inner.iter()
    }
}


#[derive(RustcDecodable, RustcEncodable)]
pub struct Record {
    pub cell: String,
    pub feature: String,
    pub expression: f64,
    pub prob: LogProb
}


pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>
}


impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Option<Self> {
        fs::File::open(path).map(|f| Reader::from_reader(f)).ok().expect("Error opening file.")
    }
}


impl<R: io::Read> Reader<R> {
    pub fn from_reader(rdr: R) -> Option<Self> {
        let mut inner = csv::Reader::from_reader(rdr).delimiter(b'\t');
        if inner.headers().unwrap() == HEADER {
            Some(Reader {
                inner: inner
            })
        }
        else {
            None
        }
    }

    pub fn cdfs(&mut self) -> CDFs {
        let mut features = collections::HashMap::new();
        for (feature, records) in self.inner.decode().map(|res| res.ok().expect("Error reading record")).group_by(
            |rec: &Record| rec.feature.clone()
        ) {
            let cdf = CDF::from_cdf(records.iter().map(|rec| {
                (rec.expression.clone(), rec.prob.clone())
            }));
            let mut cdfs = features.entry(feature).or_insert(Vec::new());
            cdfs.push(cdf);
        }
        CDFs { inner: features }
    }
}


pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>
}


impl Writer<fs::File> {
    /// Write to a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        fs::File::create(path).map(|f| Writer::from_writer(f)).unwrap()
    }
}


impl<W: io::Write> Writer<W> {
    pub fn from_writer(w: W) -> Self {
        let mut writer = Writer {
            inner: csv::Writer::from_writer(w).delimiter(b'\t')
        };
        writer.inner.write(HEADER.iter()).unwrap();

        writer
    }

    pub fn write(&mut self, cell: &str, feature: &str, cdf: &CDF) {
        for x in cdf.iter() {
            self.inner.write([
                cell,
                feature,
                &format!("{:.0}", x.0)[..],
                &format!("{}", x.1)[..]
            ].iter()).unwrap();
        }
    }
}
