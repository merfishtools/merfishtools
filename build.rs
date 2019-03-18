use std::fs;
use std::io::Write;

fn bit_count(a: usize) -> usize {
    let mut d = 0;
    let mut a = a;
    while a > 0 {
        d += a & 1;
        a >>= 1;
    }
    d
}

fn main() -> std::io::Result<()> {
    let nbits: Vec<usize> = (0..2 << (16 - 1)).map(|i| bit_count(i)).collect();
    let mut f = fs::File::create("src/model/la/hamming.rs")?;
    f.write_all(b"const _NBITS16: [usize; 65536] = [")?;
    f.write(b"\t")?;
    for (i, n) in nbits.iter().enumerate() {
        if i % 16 == 0 {
            f.write(b"\n\t")?;
        }
        f.write_fmt(format_args!("{}, ", n))?;
    }
    f.write(b"];\n")?;
    f.write(b"
pub fn hamming_distance16(a: usize, b: usize) -> usize {
    let x = a ^ b;
    _NBITS16[x]
 }")?;
    f.sync_all()?;
    Ok(())
}