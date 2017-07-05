use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


#[test]
fn test_exp() {
    assert!(Command::new("bash")
            .arg("-c")
            .arg("target/debug/merfishtools exp test/codebook/140genesData.1.txt < test/data/140genesData.1.cell0.txt")
            .spawn().unwrap().wait().unwrap().success());
}
