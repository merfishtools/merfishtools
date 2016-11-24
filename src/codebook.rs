use bit_set::BitSet;
use rand;
use rand::Rng;


pub struct Generator {
    n: u8,
    m: u8,
    dist: u8,
    is_valid: BitSet,
    visited: BitSet,
    queue: Vec<u16>
}

impl Generator {
    pub fn new(n: u8, m: u8, dist: u8) -> Self {
        assert!(n <= 16);
        assert!(dist % 2 == 0);
        Generator {
            n: n,
            m: m,
            dist: dist,
            is_valid: BitSet::with_capacity(Self::calc_max_value(n, m) as usize + 1),
            visited: BitSet::new(),
            queue: Vec::new()
        }
    }

    fn calc_max_value(n: u8, m: u8) -> u16 {
        let ones = (1 << m as u16) - 1;
        let max = ones << (n - m);
        max
        //(2usize.pow(n as u32) - 1) as u16
    }

    fn max_value(&self) -> u16 {
        Self::calc_max_value(self.n, self.m)
    }

    pub fn generate(&mut self) {
        self.is_valid.clear();
        self.queue.clear();
        let mut rng = rand::thread_rng();

        // step 1: initialize bitset and queue with one entry for each possible codeword
        for w in 0..self.max_value() + 1 {
            let valid = w.count_ones() == self.m as u32;
            if valid {
                self.queue.push(w);
                self.is_valid.insert(w as usize);
            }
        }
        println!("{:?}", self.is_valid);
        // step 2: shuffle queue
        rng.shuffle(&mut self.queue);

        let first = self.queue[self.queue.len() - 1];
        while let Some(w) = self.queue.pop() {
            // step 3: take random word
            if !self.is_valid.contains(w as usize) {
                // if word has been marked as invalid, skip it
                continue;
            }

            // step 4: mark neighborhood as invalid
            let dist = self.dist;
            println!("w={:016b}, first={:016b}", w, first);
            self.visited.clear();
            self.mark_neighborhood(w, dist - 1);
            assert!(self.is_valid.contains(first as usize));
        }
    }

    pub fn iter(&self) -> CodebookIterator {
        CodebookIterator { is_valid: &self.is_valid, w: 1 << self.m - 1, max_value: self.max_value() }
    }

    fn mark_neighborhood(&mut self, w: u16, max_dist: u8) {
        if max_dist < 2 {
            return;
        }
        self.visited.insert(w as usize);

        for i in 0..self.n {
            // swap two bits in order to keep m constant
            // find 1-bit to remove
            let ibit = 1 << i;
            if w & ibit == 0 {
                // continue if no 1-bit is removed
                continue;
            }
            for j in 0..self.n {
                // find 1-bit to add
                let jbit = 1 << j;
                if w & jbit > 0 {
                    // continue if no 1-bit would be added
                    continue;
                }
                let mut neighbor = w;
                neighbor ^= ibit;
                neighbor ^= jbit;
                if !self.visited.contains(neighbor as usize) {
                    self.is_valid.remove(neighbor as usize);
                    self.mark_neighborhood(neighbor, max_dist - 2);
                }
            }
        }
    }
}


pub struct CodebookIterator<'a> {
    is_valid: &'a BitSet,
    w: u16,
    max_value: u16
}


impl<'a> Iterator for CodebookIterator<'a> {
    type Item = u16;

    fn next(&mut self) -> Option<u16> {
        while self.w < self.max_value {
            let valid = self.is_valid.contains(self.w as usize);
            let w = self.w;
            self.w += 1;
            if valid {
                return Some(w);
            }
        }
        None
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_generate() {
        let dist = 4;
        let m = 6;
        let n = 16;
        let mut generator = Generator::new(n, m, dist);
        generator.generate();

        let words = generator.iter().collect_vec();
        println!("{}", words.len());
        assert!(words.len() > 0);
        for w in words.iter() {
            assert_eq!(w.count_ones() as u8, m);
        }
        let mut is_tight = false;
        for (a, b) in words.iter().cartesian_product(words.iter()) {
            if a == b {
                continue;
            }
            let d = (a ^ b).count_ones() as u8;
            assert!(d >= dist);
            if d == 4 {
                is_tight = true;
            }
        }
        assert!(is_tight);
    }
}
