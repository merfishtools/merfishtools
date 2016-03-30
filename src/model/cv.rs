use itertools::Itertools;

use model;


pub type CV = f64;


pub fn pmf(pmfs: &[model::expressionset::PMF]) -> model::diffexp::PMF {
    let meanvar = model::pmf::MeanVar::new(pmfs);

    let mut pmf = meanvar.iter().map(|e| model::pmf::Entry {
        value: e.standard_deviation().log2() - e.mean.log2(),
        prob: e.prob
    }).collect_vec();
    pmf.sort_by(|a, b| a.value.partial_cmp(&b.value).unwrap());

    model::diffexp::PMF::new(pmf)
}


#[cfg(test)]
mod tests {
    use super::*;

    use itertools::Itertools;
    use bio::stats::logprobs;

    use model;
    use io;

    const GENE: &'static str = "COL5A1";

    fn setup() -> Box<model::readout::Model> {
        model::readout::new_model(16, 4, 0.04, 0.1, 4, io::codebook::Reader::from_file("evaluation/codebook/140genesData.1.txt", 4).unwrap().codebook())
    }

    #[test]
    fn test_pmf() {
        let readout = setup();
        let pmfs1 = [
            model::expression::pmf(GENE, 5, 5, &readout, 100),
            model::expression::pmf(GENE, 5, 5, &readout, 100),
            model::expression::pmf(GENE, 5, 5, &readout, 100),
            model::expression::pmf(GENE, 5, 5, &readout, 100)
        ];
        let pmfs2 = [
            model::expression::pmf(GENE, 50, 50, &readout, 100),
            model::expression::pmf(GENE, 50, 50, &readout, 100),
            model::expression::pmf(GENE, 50, 50, &readout, 100),
            model::expression::pmf(GENE, 50, 50, &readout, 100)
        ];
        let pmf1 = model::expressionset::pmf(&pmfs1);
        let pmf2 = model::expressionset::pmf(&pmfs2);
        //println!("{} {}", pmfs1[0].expected_value(), pmfs2[0].expected_value());

        let pmf = pmf(&[pmf1, pmf2]);

        let total = logprobs::sum(&pmf.iter().map(|fc| fc.prob).collect_vec());

        assert!(total <= 0.0);
        assert_relative_eq!(total, 0.0, epsilon = 0.0002);
        assert_relative_eq!(2.0f64.powf(pmf.expected_value()), 1.14, epsilon = 0.01);
    }
}
