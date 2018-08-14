use std::process::Command;
use std::path::Path;

// fn test_output(result: &str, expected: &str) {
//     assert!(Command::new("cmp")
//             .arg(result)
//             .arg(expected)
//             .spawn().unwrap().wait().unwrap().success());
//     fs::remove_file(result).unwrap();
// }

fn run_cmd(cmd: &str, prefix: Option<&str>) -> bool {
    let build = if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    };

    println!("{}", cmd);

    Command::new("bash")
        .arg("-c")
        .arg(format!(
            "{prefix}target/{build}/{cmd}",
            prefix = prefix.unwrap_or(""),
            cmd = cmd,
            build = build
        ))
        .spawn()
        .unwrap()
        .wait()
        .unwrap()
        .success()
}

fn run_exp(dataset: &str, codebook: &str, params: &str) -> bool {
    let mut raw_path = format!("tests/data/{}.txt", dataset);
    if !Path::new(&raw_path).exists() {
        raw_path = format!("tests/data/{}.bin", dataset);
    }
    run_cmd(
        &format!(
            "merfishtools -v exp --seed 42 -t 1 {codebook} {raw} --estimate {est} --stats {stats} {params} > {cdf}",
            codebook=format!("tests/codebook/{}.txt", codebook),
            raw=raw_path,
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
    assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.01258132102084181 0.013360438064572837 0.019941605919281354 0.016581551788592854 0.22860571489439577 0.22269171721737882 0.2073701214077358 0.02231576033300221 0.017689514891827997 0.02054328192011981 0.02101414200654054 0.03938992313394839 0.023995219270935196 0.03557793748730565 0.17007370565085334 0.03569158826173908 --p1 0.5912118122316706 0.5063863207089024 0.42338772914919315 0.43282699715628403 0.17843503208026912 0.18766044306365068 0.1989449380102563 0.42262699427628003 0.4930707032134955 0.4034046655581933 0.4610572590564446 0.30553305919201673 0.2931349183480258 0.2634814096058736 0.16619112507647252 0.36410395236487053"));
    //assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.016 0.015 0.021 0.015 0.051 0.065 0.062 0.026 0.019 0.021 0.023 0.034 0.02 0.029 0.07 0.027 --p1 0.05 0.04 0.065 0.052 0.177 0.135 0.155 0.07 0.053 0.057 0.066 0.09 0.069 0.07 0.145 0.103"));
    //assert!(run_exp("140genesData.1.cell0", "140genesData.1", "--p0 0.1 --p1 0.6"));
}

#[test]
fn test_exp_mhd4_thbs_bias() {
    assert!(run_exp("140genesData.1.cell24", "140genesData.1", "--p0 0.01258132102084181 0.013360438064572837 0.019941605919281354 0.016581551788592854 0.22860571489439577 0.22269171721737882 0.2073701214077358 0.02231576033300221 0.017689514891827997 0.02054328192011981 0.02101414200654054 0.03938992313394839 0.023995219270935196 0.03557793748730565 0.17007370565085334 0.03569158826173908 --p1 0.5912118122316706 0.5063863207089024 0.42338772914919315 0.43282699715628403 0.17843503208026912 0.18766044306365068 0.1989449380102563 0.42262699427628003 0.4930707032134955 0.4034046655581933 0.4610572590564446 0.30553305919201673 0.2931349183480258 0.2634814096058736 0.16619112507647252 0.36410395236487053"));
}

#[test]
fn test_exp_mhd2() {
    assert!(run_exp("1001genesData.3.cell0", "1001genesData", ""));
}

/// This is a longrunning test that should be activated only on purpose.
#[test]
#[ignore]
fn test_exp_mhd2_8() {
    assert!(run_exp(
        "simulated-MHD2-8.25.all",
        "simulated-MHD2-8",
        "--p0 0.005 --p1 0.01"
    ));
}

#[test]
fn test_exp_mhd4_sim() {
    assert!(run_exp(
        "simulated-MHD4.35.1",
        "simulated-MHD4",
        "--p0 0.04 --p1 0.1"
    ));
}

#[test]
fn test_exp_mhd4v2_real() {
    assert!(run_exp(
        "real-MHD4v2.small-example",
        "real-MHD4v2",
        "--p0 0.04 --p1 0.1"
    ));
}

#[test]
fn test_estimate_error_rates_real() {
    assert!(run_cmd(
        "merfishtools -v est-error-rates tests/codebook/140genesData.1.txt > tests/results/140genesData.error-rates.tsv",
        Some("grep -P '^3\\t' tests/data/140genesData.readouts.txt | cut -f2,3,4 | ")
    ));
}

#[test]
fn test_estimate_error_rates_simulated() {
    assert!(run_cmd(
        "merfishtools -v est-error-rates tests/codebook/140genesData.1.txt < tests/data/simulated-MHD4.35.readouts.txt > tests/results/simulated-MHD4.35.error-rates.tsv",
        None
    ));
}
