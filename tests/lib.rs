use std::process::Command;
use std::fs;


fn test_output(result: &str, expected: &str) {
    assert!(Command::new("cmp")
            .arg(result)
            .arg(expected)
            .spawn().unwrap().wait().unwrap().success());
    fs::remove_file(result).unwrap();
}


fn run_cmd(cmd: &str, prefix: Option<&str>) -> bool {
    let build = if cfg!(debug_assertions) { "debug" } else { "release" };

    println!("{}", cmd);

    Command::new("bash")
            .arg("-c")
            .arg(format!(
                "{prefix}target/{build}/{cmd}",
                prefix=prefix.unwrap_or(""),
                cmd=cmd,
                build=build
            ))
            .spawn().unwrap().wait().unwrap().success()
}


fn run_exp(dataset: &str, codebook: &str, params: &str) -> bool {
    run_cmd(
        &format!(
            "merfishtools -v exp -t 1 {codebook} --estimate {est} --stats {stats} {params} < {raw} > {cdf}",
            codebook=format!("tests/codebook/{}.txt", codebook),
            raw=format!("tests/data/{}.txt", dataset),
            cdf=format!("tests/results/{}.txt", dataset),
            est=format!("tests/results/{}.est.txt", dataset),
            stats=format!("tests/results/{}.stats.txt", dataset),
            params=params
        ),
        None
    )
}


#[test]
fn test_exp_mhd4_cell0() {
    assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.016 0.015 0.021 0.015 0.051 0.065 0.062 0.026 0.019 0.021 0.023 0.034 0.02 0.029 0.07 0.027 --p1 0.05 0.04 0.065 0.052 0.177 0.135 0.155 0.07 0.053 0.057 0.066 0.09 0.069 0.07 0.145 0.103"));
    //assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.1 --p1 0.6"));
}


#[test]
fn test_exp_mhd4_thbs_bias() {
    assert!(run_exp("140genesData.1.cell24", "140genesData.1", "--p0 0.016 0.015 0.021 0.015 0.051 0.065 0.062 0.026 0.019 0.021 0.023 0.034 0.02 0.029 0.07 0.027 --p1 0.05 0.04 0.065 0.052 0.177 0.135 0.155 0.07 0.053 0.057 0.066 0.09 0.069 0.07 0.145 0.103"));
}


#[test]
fn test_exp_mhd2() {
    assert!(run_exp("1001genesData.3.cell0", "1001genesData", ""));
}


/// This is a longrunning test that should be activated only on purpose.
fn test_exp_mhd2_8() {
    assert!(run_exp("simulated-MHD2-8.25.all", "simulated-MHD2-8", "--p0 0.005 --p1 0.01"));
}

#[test]
fn test_exp_mhd4_sim() {
    assert!(run_exp("simulated-MHD4.35.1", "simulated-MHD4", "--p0 0.04 --p1 0.1"));
}

#[test]
fn test_estimate_error_rates_real() {
    assert!(run_cmd(
        "merfishtools -v est-error-rates tests/codebook/140genesData.1.txt > tests/results/140genesData.error-rates.tsv",
        Some("grep -P '^1\\t' tests/data/140genesData.readouts.txt | ")
    ));
}


#[test]
fn test_estimate_error_rates_simulated() {
    assert!(run_cmd(
        "merfishtools -v est-error-rates tests/codebook/140genesData.1.txt < tests/data/simulated-MHD4.35.readouts.txt > tests/results/simulated-MHD4.35.error-rates.tsv",
        None
    ));
}
