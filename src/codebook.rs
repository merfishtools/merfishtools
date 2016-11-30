
use std::collections::VecDeque;
use bit_set::BitSet;
use rand;
use rand::Rng;

use itertools::Itertools;
use ndarray::prelude::*;


pub fn count_ones<'a, I: Iterator<Item=&'a u8>>(bits: I) -> u8 {
    bits.fold(0, |count, &b| count + b)
}


/// Generate MHD4 code with m 1-bits.
pub fn generate_mhd4(m: u8) -> Vec<Vec<u8>> {
    let gen_mat = Array::from_shape_vec((11, 16), vec![
        1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,
        1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,
        0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,0,
        1,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,
        1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,
        0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,
        0,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,
        1,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,
        0,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0,
        1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1
    ]).unwrap();

    let mut words = Vec::new();
    for i in 0..2u32.pow(11) - 1 {
        let mut vec = Vec::new();
        for j in 0..11 {
            vec.push(((i >> j) & 1) as u8);
        }
        let vec = Array::from_shape_vec((1, 11), vec).unwrap();
        let word = vec.dot(&gen_mat) % 2;

        if count_ones(word.iter()) == m {
            assert_eq!(word.shape(), &[1, 16]);
            words.push(word.into_raw_vec());
            //println!("{}", word.iter().map(|b| format!("{}", b)).join(""));
        }
    }

    words
}


/// Generate MHD2 code with n bits in total and m 1-bits.
pub fn generate_mhd2(n: u8, m: u8) -> Vec<Vec<u8>> {
    (0..n as usize).combinations(m as usize).map(|bits| {
        let mut word = vec![0; n as usize];
        for i in bits {
            word[i] = 1;
        }
        word
    }).collect_vec()
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    fn hamming_dist<'a, I: Iterator<Item=&'a u8>>(a: I, b: I) -> u8 {
        a.zip(b).fold(0, |d, (&a, &b)| d + (if a != b {1} else {0}))
    }


    fn test_words(words: &[Vec<u8>], n: u8, m: u8, dist: u8) {
        for w in words.iter() {
            assert_eq!(count_ones(w.iter()) as u8, m);
            assert_eq!(w.len(), n as usize);

        }

        let mut is_tight = false;
        for (a, b) in words.iter().cartesian_product(words.iter()) {
            if a == b {
                continue;
            }

            let d = hamming_dist(a.iter(), b.iter());
            assert!(d >= dist, format!("{} < {}", d, dist));
            if d == dist {
                is_tight = true;
            }
        }
        assert!(is_tight);
    }


    #[test]
    fn test_generate_mhd2() {
        let dist = 2;
        let m = 4;
        let n = 14;
        let words = generate_mhd2(n, m);
        assert_eq!(words.len(), 1001);
        test_words(&words, n, m, dist);
    }


    #[test]
    fn test_generate_mhd4() {
        let dist = 4;
        let m = 6;
        let n = 16;
        let words = generate_mhd4(m);
        println!("{}", words.len());
        assert_eq!(words.len(), 448);
        test_words(&words, n, m, dist);
    }
}
