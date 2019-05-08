use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use derive_new::new;
use failure::Error;
use itertools::*;
use itertools::Itertools;
use ndarray::{Array1, Axis};
use ndarray::s;
use serde::{Deserialize, Serialize};
use structopt::StructOpt;

use merfishtools::model::la::expression::{estimate_errors, estimate_expression};
use merfishtools::model::la::hamming::_NBITS16;
use merfishtools::model::la::matrix::csr_error_matrix;
use merfishtools::simulation::*;
use serde::de::Unexpected::StructVariant;
use std::cmp::Ordering;

mod binary {
    use serde::{Deserializer, Serializer};
    use serde::de::Deserialize;

    use crate::Barcode;

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
            acc | ((c.to_digit(2).unwrap() as u16) << i)
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

type Expr = Array1<f32>;
type Barcode = u16;

fn _read_records<R: std::io::Read, T: for<'de> Deserialize<'de>>(
    reader: R,
) -> Result<Vec<T>, Error> {
    let mut raw_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(reader);
    Ok(raw_reader.deserialize::<T>().map(|r| r.unwrap()).collect())
}

#[derive(StructOpt)]
struct SimulationParams {
    raw_counts: String,
    obs_counts: String,
}

fn main() -> Result<(), Error> {
    let args: SimulationParams = StructOpt::from_args();
    let raw_path = File::open(args.raw_counts)?;
    let obs_path = File::open(args.obs_counts)?;
    let raw_records = _read_records::<_, RecordRaw>(raw_path)?;
    let raw_counts: HashMap<usize, HashMap<Barcode, usize>> = raw_records
        .iter()
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
        .collect();
    let barcodes: std::collections::HashSet<Barcode> = raw_counts.iter().flat_map(|(_cell, counts)| counts.keys().cloned()).collect();
    let obs_records = _read_records::<_, RecordObserved>(obs_path)?;
    let obs_derived: HashMap<usize, HashMap<Barcode, Vec<(Barcode, usize)>>> = obs_records
        .into_iter()
        .group_by(|r| r.cell)
        .into_iter()
        .map(|(cell, recs)| (cell, recs.into_iter().collect_vec()))
        .map(|(cell, recs)| {
            let foo: HashMap<Barcode, Vec<(Barcode, usize)>> = recs
                .iter()
                .group_by(|&r| r.readout)
                .into_iter()
                .map(|(readout, subrecs)| {
                    (
                        readout,
                        subrecs
                            .map(|subrec| (subrec.barcode, subrec.count))
                            .collect_vec(),
                    )
                })
                .collect();
            (cell, foo)
        })
        .collect();
    let obs_counts: HashMap<usize, HashMap<Barcode, usize>> = obs_derived
        .iter()
        .map(|(&cell, readouts)| {
            (
                cell,
                readouts
                    .iter()
                    .map(|(&readout, barcodes)| (readout, barcodes.iter().map(|b| b.1).sum()))
                    .collect(),
            )
        })
        .collect();
    const NUM_BITS: usize = 16;
    const NUM_CODES: usize = 1 << NUM_BITS;
    let mut e = [[0.; NUM_BITS], [0.; NUM_BITS]];
    let p0 = vec![0.05; NUM_BITS];
    let p1 = vec![0.09; NUM_BITS];
    e[0].iter_mut()
        .zip(p0.iter().cloned().map(f64::from))
        .for_each(|(v, p)| *v = p as f32);
    e[1].iter_mut()
        .zip(p1.iter().cloned().map(f64::from))
        .for_each(|(v, p)| *v = p as f32);

    let mut x_est: Expr = Array1::<f32>::zeros(NUM_CODES);
    let mut y: Expr = Array1::<f32>::zeros(NUM_CODES);
    for (_cell, feature_counts) in &obs_counts {
        for (&barcode, &count) in feature_counts {
            let bc = barcode as usize;
            if _NBITS16[bc] == 4 {
                x_est[bc] += count as f32;
            } else if _NBITS16[bc] == 3 || _NBITS16[bc] == 5 {
                for i in 0..NUM_BITS {
                    let modified = barcode ^ (1 << i);
                    if _NBITS16[modified as usize] == 4 && raw_counts.get(&_cell).unwrap().contains_key(&modified) {
                        x_est[modified as usize] += count as f32;
                    }
                }
            }
        }
    }
    let x_est_initial = x_est.clone();
    for (_cell, feature_counts) in &obs_counts {
        for (&barcode, &count) in feature_counts {
            y[barcode as usize] += count as f32;
        }
    }
    let yv = y.view();
    let mut x: Expr = x_est;
//    dbg!(&x);
    let xx = raw_counts
        .iter()
        .flat_map(|(_cell, counts)| counts.into_iter().collect::<Vec<(&Barcode, &usize)>>())
        .map(|(&a, &b)| (a, b)).collect_vec();
    let mut x_true: Expr = Expr::zeros(NUM_CODES);
    for (bc, c) in xx {
        x_true[bc as usize] += c as f32;
    }
    fn barcode_mask(a: &Expr, barcodes: &std::collections::HashSet<Barcode>) -> Expr {
        let mut e = Expr::zeros(NUM_CODES);
        for &bc in barcodes {
            e[bc as usize] = a[bc as usize];
        }
        e
    }
//    for hamming_dist in 1..=4 {
//        println!("\n=== {:?} ===", hamming_dist);
//        e = estimate_errors(x.view(), yv, &e, hamming_dist).expect("Failed estimating errors");
//        dbg!(&e);
//        let mat = csr_error_matrix(&e, hamming_dist.max(2));
//        x = estimate_expression(&mat, x.view(), yv, hamming_dist, false).expect("Failed estimating expression");
//        dbg!(mse(&x_true, &x));
////        x.iter_mut().enumerate().filter(|(i, _)| !barcodes.contains(&i)).for_each(|(_, v)| *v /= 1.1);
//        x.iter_mut().filter(|&&mut v| v <= 1.).for_each(|v| *v = 0.);
//        dbg!(mse(&x_true, &x));
////        dbg!(mse(&x_true, &x));
////        dbg!(&e);
//
//    }
    for hamming_dist in 2..=4 {
        let mat = csr_error_matrix(&e, hamming_dist.max(2));
        println!("Estimating x …");
        x = estimate_expression(&mat, x.view(), yv, hamming_dist, false).expect("Failed estimating expression");
        dbg!(mse(&x_true, &x));
//        if hamming_dist == 1 {
            println!("Estimating E …");
            e = estimate_errors(x.view(), yv, &e, hamming_dist).expect("Failed estimating errors");
            dbg!(&e);
//        }
    }
    dbg!(mse(&x_true, &x));
    dbg!(mse(&x_true, &barcode_mask(&x, &barcodes)));

    dbg!(mse(&x_true, &x_est_initial));
    dbg!(mse(&x_true, &barcode_mask(&x_est_initial, &barcodes)));

    dbg!(mse(&x_true, &y));
    dbg!(mse(&x_true, &barcode_mask(&y, &barcodes)));

    Ok(())
}

fn mse(a: &Expr, b: &Expr) -> f32 {
    let diff: Expr = a - b;
    diff.mapv(|x| x.powi(2)).sum() / (a.len() as f32)
}
