use std::fs;
use std::io::Write;

fn main() -> std::io::Result<()> {
    if fs::File::open("src/model/la/mask.rs").is_ok() {
        return Ok(());
    }
    let mut f = fs::File::create("src/model/la/mask.rs")?;
    f.write_all(b"pub const _DIFFZERO_MASK_16: [[u32; 65536]; 16] = [\n")?;
    for nbits in 0..=16 {
        let nonzero_mask: Vec<u32> = (0..1 << 16).map(|v: u32| if v.count_ones() <= nbits { v } else { 0 }).collect();
        f.write(b"    \n[")?;
        for (i, v) in nonzero_mask.iter().enumerate() {
            if i % 16 == 0 {
                f.write(b"\n     ")?;
            }
            f.write_fmt(format_args!("{}, ", v))?;
        }
        f.write(b"],\n")?;
    }
    f.write(b"];\n")?;
    f.sync_all()?;
    Ok(())
}