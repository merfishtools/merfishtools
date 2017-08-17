use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


fn run_exp(dataset: &str, codebook: &str, params: &str) -> bool {
    Command::new("bash")
            .arg("-c")
            .arg(format!(
                "target/debug/merfishtools -v exp -t 1 {codebook} --estimate {est} {params} < {raw} > {cdf}",
                codebook=format!("tests/codebook/{}.txt", codebook),
                raw=format!("tests/data/{}.txt", dataset),
                cdf=format!("tests/results/{}.txt", dataset),
                est=format!("tests/results/{}.est.txt", dataset),
                params=params
            ))
            .spawn().unwrap().wait().unwrap().success()
}


#[test]
fn test_exp_mhd4_cell0() {
    assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.016 0.015 0.021 0.015 0.051 0.065 0.062 0.026 0.019 0.021 0.023 0.034 0.02 0.029 0.07 0.027 --p1 0.05 0.04 0.065 0.052 0.177 0.135 0.155 0.07 0.053 0.057 0.066 0.09 0.069 0.07 0.145 0.103"));
}


#[test]
fn test_exp_mhd4_thbs_bias() {
    assert!(run_exp("140genesData.1.cell24", "140genesData.1", "--p0 0.016 0.015 0.021 0.015 0.051 0.065 0.062 0.026 0.019 0.021 0.023 0.034 0.02 0.029 0.07 0.027 --p1 0.05 0.04 0.065 0.052 0.177 0.135 0.155 0.07 0.053 0.057 0.066 0.09 0.069 0.07 0.145 0.103"));
}


#[test]
fn test_exp_mhd2() {
    assert!(run_exp("1001genesData.3.cell0", "1001genesData", ""));
}
