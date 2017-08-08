use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


fn run_exp(inpath: &str, codebook: &str) -> bool {
    Command::new("bash")
            .arg("-c")
            .arg(format!(
                "target/release/merfishtools -v exp {codebook} < {raw}",
                codebook=codebook,
                raw=inpath
            ))
            .spawn().unwrap().wait().unwrap().success()
}


#[test]
fn test_exp_mhd4() {
    assert!(run_exp("tests/data/140genesData.1.cell0.txt", "tests/codebook/140genesData.1.txt"));
}


#[test]
fn test_exp_mhd2() {
    assert!(run_exp("tests/data/1001genesData.3.cell0.txt", "tests/codebook/1001genesData.txt"));
}
