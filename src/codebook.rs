// Copyright 2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bit_vec::BitVec;
use itertools::Itertools;
use ndarray::prelude::*;

pub type Word = BitVec;

/// Generate MHD4 code with m 1-bits.
pub fn generate_mhd4(m: u8) -> Vec<Word> {
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

    let mut words = Vec::new();
    for i in 0..2u32.pow(11) - 1 {
        let mut vec = Vec::new();
        for j in 0..11 {
            vec.push(((i >> j) & 1) as u8);
        }
        let vec = Array::from_shape_vec((1, 11), vec).unwrap();
        let word = vec.dot(&gen_mat) % 2;
        let word = BitVec::from_fn(16, |i| word[[0, i]] == 1);

        if word.iter().filter(|x| *x).count() == m as usize {
            words.push(word);
        }
    }

    words
}

/// Generate MHD2 code with n bits in total and m 1-bits.
pub fn generate_mhd2(n: u8, m: u8) -> Vec<Word> {
    (0..n as usize)
        .combinations(m as usize)
        .map(|bits| {
            let mut word = BitVec::from_elem(n as usize, false);
            for i in bits {
                word.set(i, true);
            }
            word
        })
        .collect_vec()
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use super::*;

    fn hamming_dist(a: &Word, b: &Word) -> u8 {
        assert_eq!(a.len(), b.len());
        (0..a.len()).fold(0, |d, i| {
            d + (if a.get(i).unwrap() != b.get(i).unwrap() {
                1
            } else {
                0
            })
        })
    }

    fn test_words(words: &[Word], dist: u8) {
        let mut is_tight = false;
        for (a, b) in words.iter().cartesian_product(words.iter()) {
            if a == b {
                continue;
            }

            let d = hamming_dist(&a, &b);
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
        test_words(&words, dist);
    }

    #[test]
    fn test_generate_mhd4() {
        let dist = 4;
        let m = 6;
        let words = generate_mhd4(m);
        println!("{}", words.len());
        assert_eq!(words.len(), 448);
        test_words(&words, dist);
    }
}
