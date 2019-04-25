use core::fmt::Pointer;
use std::collections::HashMap;

use bio::stats::Prob;
use derive_new::*;
use failure::Error;
use itertools::Itertools;
use maplit::hashmap;
use ndarray::Array;
use rand::prelude::*;
use rand::{Rng, SeedableRng};
use serde::{Serialize, Serializer};

use crate::model::la::hamming::_NBITS16;
use std::path::Path;

type Barcode = u16;

#[rustfmt::skip]
pub fn generate_barcodes(bits: usize, hamming_distance: usize) -> Vec<Barcode> {
    // MHD4
    if bits == 16 && hamming_distance == 4 {
        let gen_mat = Array::from_shape_vec(
            (11, 16),
            vec![
                1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            ],
        ).unwrap();

        let barcodes: Vec<Barcode> = (0..2u32.pow(11) - 1)
            .map(|i| {
                let vec = (0..11).map(|j| ((i >> j) & 1) as u8).collect();
                let vec = Array::from_shape_vec((1, 11), vec).unwrap();
                let word = vec.dot(&gen_mat) % 2;
                let word: u16 = word.iter().fold(0, |acc, v| (acc << 1) | (*v as u16));
                word
            })
            .collect();
        return barcodes;
    }
    unimplemented!()
}

pub fn generate_raw_counts(
    barcodes: &[Barcode],
    lambda: f64,
    rng: &mut StdRng,
) -> HashMap<Barcode, usize> {
    let num_codes = barcodes.len();
    let poisson = rand::distributions::Poisson::new(lambda);

    poisson
        .sample_iter(rng)
        .enumerate()
        .take(num_codes)
        .map(|(i, v)| (barcodes[i], v as usize))
        .collect()
}

pub fn generate_erroneous_counts(
    raw_counts: &HashMap<Barcode, usize>,
    rng: &mut StdRng,
    p0: &[f32],
    p1: &[f32],
    max_hamming_distance: usize,
) -> HashMap<Barcode, HashMap<Barcode, usize>> {
    assert!(p0.len() <= 16);
    assert_eq!(p0.len(), p1.len());
    let num_bits = p0.len();
    let uniform = rand::distributions::Uniform::new(0f32, 1f32);
    let mut derived_counts: HashMap<Barcode, HashMap<Barcode, usize>> = raw_counts
        .to_owned()
        .iter()
        .map(|(&k, &v)| (k, hashmap! {0 => v}))
        .collect();

    raw_counts.iter().for_each(|(barcode, &count)| {
        (0..count).for_each(|_c| {
            let mut b = *barcode;
            let mut flip = 0u16;
            // choose max_hamming_distance different positions to flip
            (0..num_bits)
                .choose_multiple(rng, max_hamming_distance)
                .iter()
                .for_each(|&idx| {
                    let bit = (barcode >> idx) & 1;
                    let p = match bit {
                        0 => p0[idx],
                        1 => p1[idx],
                        _ => panic!("bit can only be 0 or 1"),
                    };
                    // flip (or don't) according to error rates
                    if rng.sample(uniform) < p {
                        flip ^= 1 << idx;
                        b ^= 1 << idx;
                    }
                });
            if flip != 0 {
                let c = derived_counts.entry(b).or_insert(hashmap! {flip => 0});
                *c.entry(flip).or_insert(0) += 1;
                derived_counts.entry(*barcode).and_modify(|v| {
                    v.entry(0).and_modify(|x| *x -= 1);
                });
            }
        });
    });
    derived_counts
}

mod binary {
    use serde::de::Deserialize;
    use serde::{Deserializer, Serializer};

    use crate::simulation::Barcode;

    pub(crate) fn serialize<S>(barcode: &Barcode, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let repr = format!("{:016b}", barcode);
        serializer.serialize_str(&repr)
    }

    pub(crate) fn deserialize<'de, D>(deserializer: D) -> Result<Barcode, D::Error>
    where
        D: Deserializer<'de>,
    {
        let repr: String = String::deserialize(deserializer)?;
        let barcode = repr.chars().rev().enumerate().fold(0u16, |acc, (i, c)| {
            acc + ((c.to_digit(2).unwrap() as u16) << i)
        });
        Ok(barcode)
    }
}

#[derive(Serialize, Deserialize, new)]
struct RecordRaw {
    cell: usize,
    #[serde(with = "binary")]
    barcode: Barcode,
    count: usize,
}

#[derive(Serialize, Deserialize, new)]
struct RecordObserved {
    cell: usize,
    #[serde(with = "binary")]
    readout: Barcode,
    count: usize,
    #[serde(with = "binary")]
    barcode: Barcode,
    errors: usize,
}

fn _read_raw_counts<R: std::io::Read>(
    reader: R,
) -> Result<HashMap<usize, HashMap<Barcode, usize>>, Error> {
    let mut raw_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    Ok(raw_reader
        .deserialize::<RecordRaw>()
        .map(|r| r.unwrap())
        .group_by(|r| r.cell)
        .into_iter()
        .map(|(cell, records)| {
            (
                cell,
                records
                    .into_iter()
                    .map(|r| (r.barcode, r.count))
                    .collect::<HashMap<Barcode, usize>>(),
            )
        })
        .collect())
}

fn _generate_raw_counts(
    num_cells: usize,
    barcodes: &[Barcode],
    lambda: f64,
    rng: &mut StdRng,
) -> HashMap<usize, HashMap<Barcode, usize>> {
    (0..num_cells)
        .map(|cell| (cell, generate_raw_counts(&barcodes, lambda, rng)))
        .collect()
}

fn _write_raw_counts<S: std::io::Write>(
    writer: &mut csv::Writer<S>,
    raw_counts: &HashMap<usize, HashMap<Barcode, usize>>,
) -> Result<(), Error> {
    for (cell, counts) in raw_counts {
        for (&barcode, &count) in counts {
            let record = RecordRaw::new(*cell, barcode, count);
            writer.serialize(record)?;
        }
        writer.flush()?;
    }
    Ok(())
}

pub fn simulate_raw_counts(
    bits: u8,
    min_hamming_distance: u8,
    raw_expression_path: Option<String>,
    set_bits: Option<u8>,
    num_cells: Option<usize>,
    num_barcodes: Option<usize>,
    lambda: f64,
    seed: Option<u64>,
) -> Result<(), Error> {
    let num_bits = bits as usize;
    let num_cells = num_cells.unwrap_or(std::usize::MAX);
    let mut rng = if let Some(seed) = seed {
        StdRng::seed_from_u64(seed)
    } else {
        StdRng::from_entropy()
    };

    let mut barcodes: Vec<Barcode> = generate_barcodes(num_bits, min_hamming_distance as usize);
    barcodes = if let Some(s) = set_bits {
        // TODO can be replaced with `drain_filter` in future rust releases
        barcodes
            .iter()
            .cloned()
            .filter(|b| _NBITS16[*b as usize] == s as usize)
            .collect()
    } else {
        barcodes
    };
    if let Some(num_barcodes) = num_barcodes {
        if num_barcodes > barcodes.len() {
            // TODO use failure and abort instead?
            warn!(
                "Desired number of barcodes ({}) > number of available barcodes ({}).",
                num_barcodes,
                barcodes.len()
            );
        }
        barcodes.shuffle(&mut rng);
        barcodes.truncate(num_barcodes);
    }

    let outlet = _writer(raw_expression_path);

    let raw_counts = _generate_raw_counts(num_cells, &barcodes, lambda, &mut rng);
    let mut raw_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(outlet);
    _write_raw_counts(&mut raw_writer, &raw_counts)?;
    Ok(())
}

fn _writer(sink: Option<String>) -> impl std::io::Write {
    let outlet: Box<std::io::Write> = match sink {
        None => Box::new(std::io::stdout()),
        Some(sink_path) => {
            let sink_path = Path::new(&sink_path);
            if sink_path.exists() {
                warn!("{:?} already exists, exiting.", sink_path);
                std::process::exit(1);
            }
            Box::new(std::fs::File::create(sink_path).unwrap())
        }
    };
    outlet
}

fn _reader(source: Option<String>) -> impl std::io::Read {
    let inlet: Box<std::io::Read> = match source {
        None => Box::new(std::io::stdin()),
        Some(source_path) => {
            let source_path = Path::new(&source_path);
            if !source_path.exists() {
                panic!("File does not exist {:?}", source_path);
            }
            Box::new(std::fs::File::open(source_path).unwrap())
        }
    };
    inlet
}

pub fn simulate_observed_counts(
    raw_expression_path: Option<String>,
    ecc_expression_path: Option<String>,
    p0: Vec<f64>,
    p1: Vec<f64>,
    seed: Option<u64>,
) -> Result<(), Error> {
    let p0: Vec<f32> = p0.iter().map(|&v| v as f32).collect();
    let p1: Vec<f32> = p1.iter().map(|&v| v as f32).collect();

    let mut rng = if let Some(seed) = seed {
        StdRng::seed_from_u64(seed)
    } else {
        StdRng::from_entropy()
    };

    let ecc_out = _writer(ecc_expression_path);
    let mut ecc_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(ecc_out);

    let raw_counts = _read_raw_counts(_reader(raw_expression_path))?;
    for (&cell, records) in &raw_counts {
        let derived_counts = generate_erroneous_counts(records, &mut rng, &p0, &p1, p0.len());
        for (barcode, errcount) in derived_counts.into_iter().sorted_by_key(|v| v.0) {
            for (flip, count) in errcount {
                let original_barcode = barcode ^ flip;
                let errs = _NBITS16[flip as usize];
                let record = RecordObserved::new(cell, barcode, count, original_barcode, errs);
                ecc_writer.serialize(record)?;
            }
        }
        ecc_writer.flush()?;
    }
    ecc_writer.flush()?;
    Ok(())
}
