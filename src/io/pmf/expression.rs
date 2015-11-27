use std::io;
use std::fs;
use std::path::Path;
use std::collections;

use itertools::Itertools;
use csv;

use bio::stats::logprobs::LogProb;

use model::expression::PMF;


pub struct PMFs {
    inner: collections::HashMap<String, Vec<PMF>>
}


impl PMFs {
    pub fn contains_feature(&self, feature: &str) -> bool {
        self.inner.contains_key(feature)
    }

    pub fn get(&self, feature: &str) -> Option<&Vec<PMF>> {
        self.inner.get(feature)
    }

    pub fn features(&self) -> collections::hash_map::Keys<String, Vec<PMF>>  {
        self.inner.keys()
    }

    pub fn iter(&self) -> collections::hash_map::Iter<String, Vec<PMF>> {
        self.inner.iter()
    }
}


#[derive(RustcDecodable, RustcEncodable)]
pub struct Record {
    pub experiment: String,
    pub cell: String,
    pub feature: String,
    pub expression: u32,
    pub prob: LogProb
}


pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>
}


impl Reader<fs::File> {
    /// Read from a given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Self {
        fs::File::open(path).map(|f| Reader::from_reader(f)).ok().expect("Error opening file.")
    }
}


impl<R: io::Read> Reader<R> {
    pub fn from_reader(rdr: R) -> Self {
        Reader {
            inner: csv::Reader::from_reader(rdr).delimiter(b'\t')
        }
    }

    pub fn pmfs(&mut self) -> PMFs {
        let mut features = collections::HashMap::new();
        for (feature, records) in self.inner.decode().map(|res| res.ok().expect("Error reading record")).group_by(
            |rec: &Record| rec.feature.clone()
        ) {
            let pmf = PMF::new(records.iter().map(|rec| (rec.expression.clone(), rec.prob.clone())).collect_vec());
            if !features.contains_key(&feature) {
                features.insert(feature.clone(), Vec::new());
            }
            let pmfs = features.get_mut(&feature).unwrap();
            pmfs.push(pmf);
        }
        PMFs { inner: features }
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
        writer.inner.write(["expmnt", "cell", "feat", "expr", "prob"].iter()).unwrap();

        writer
    }

    pub fn write(&mut self, experiment: &str, cell: &str, feature: &str, pmf: &PMF) {
        let mut record = Record {
            experiment: experiment.to_owned(),
            cell: cell.to_owned(),
            feature: feature.to_owned(),
            expression: 0,
            prob: 0.0
        };

        for &(x, prob) in pmf.iter() {
            record.expression = x;
            record.prob = prob;
            self.inner.encode(&record).unwrap();
        }
    }
}
