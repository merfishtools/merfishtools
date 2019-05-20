pub fn hamming_distance16(a: usize, b: usize) -> usize {
    let x = a ^ b;
    x.count_ones() as usize
}