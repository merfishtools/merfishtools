use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


fn run_exp(dataset: &str, codebook: &str) -> bool {
    Command::new("bash")
            .arg("-c")
            .arg(format!(
                "target/release/merfishtools -v exp {codebook} --estimate {est} < {raw} > {cdf}",
                codebook=format!("tests/codebook/{}.txt", codebook),
                raw=format!("tests/data/{}.txt", dataset),
                cdf=format!("tests/results/{}.txt", dataset),
                est=format!("tests/results/{}.est.txt", dataset)
            ))
            .spawn().unwrap().wait().unwrap().success()
}


#[test]
fn test_exp_mhd4() {
    assert!(run_exp("140genesData.1.cell0", "140genesData.1"));
}


#[test]
fn test_exp_mhd2() {
    assert!(run_exp("1001genesData.3.cell0", "1001genesData"));
}
