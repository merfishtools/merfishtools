import pandas as pd
sys.path.insert(0, "../")
import merfishtools as mt


configfile: "config.yaml"


merfishtools = "../target/release/merfishtools"


contexts = ["paper"]
datasets = ["140genesData"]


def experiments(dataset):
    return range(1, config["datasets"][dataset]["experiments"] + 1)


def matrices(dataset, type="expressions", settings="default"):
    suffix = "." if type == "counts" else ".{}.".format(settings)
    return expand("{type}/{dataset}.{experiment}.all{suffix}matrix.txt",
                  dataset=dataset,
                  type=type,
                  suffix=suffix,
                  experiment=experiments(dataset))


rule all:
    input:
        expand([
            "results/{context}/foldchange_cdf/{dataset}.1.A-vs-B.{gene}.default.foldchange_cdf.pdf",
            "results/{context}/expression_pmf/{dataset}.1.cell0.{gene}.default.expression_pmf.pdf"
        ], gene=config["genes"], context=contexts, dataset=datasets),
        expand([
            "results/{context}/{dataset}.default.expression_dist.pdf",
            "results/{context}/{dataset}.default.overdispersion.pdf",
            "results/{context}/{dataset}.default.correlation.pdf"
        ], context=contexts, dataset=datasets),
        expand(["results/{context}/{dataset}.{type}.default.pca.pdf",
                "results/{context}/{dataset}.{type}.default.qqplot.pdf"],
               context=contexts, dataset=datasets,
               type=["expressions", "normalized_expressions"]),
        expand("results/{context}/simulation-MHD{dist}/MHD{dist}.{plot}.default.pdf", plot=["scatter", "error"], context=contexts, dist=[2, 4]),
        expand("results/{context}/default.dataset_correlation.pdf", context=contexts)


rule format:
    input:
        "data/{dataset}.csv.bz2"
    output:
        "data/{dataset}.{experiment}.{group}.txt"
    script:
        "scripts/format-dataset.py"



rule raw_counts:
    input:
        "data/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.txt"
    script:
        "scripts/raw-counts.py"


rule expressions:
    input:
        "data/{dataset}.{experiment}.{group}.txt"
    output:
        pmf="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        est="expressions/{dataset}.{experiment}.{group}.{settings}.est.txt",
    params:
        dist=lambda wildcards: config["datasets"][wildcards.dataset]["dist"],
        bits=lambda wildcards: config["datasets"][wildcards.dataset]["N"]
    benchmark:
        "bench/exp/{dataset}.{settings}.txt"
    threads: 8
    shell:
        "{merfishtools} exp --hamming-dist {params.dist} -N {params.bits} "
        "--estimate {output.est} -t {threads} < {input} > {output.pmf}"


rule normalize_pmf:
    input:
        pmf="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings,(default)}.txt"
    run:
        pmf = pd.read_table(input.pmf, index_col=0)
        scales = pd.read_table(input.scales, index_col=0, squeeze=True, header=None)
        pmf["expr"] *= scales[int(wildcards.experiment)]
        pmf.to_csv(output[0], sep="\t")


rule diffexp:
    input:
        "expressions/{dataset}.{experiment1}.{group1}.{settings}.txt",
        "expressions/{dataset}.{experiment2}.{group2}.{settings}.txt"
    output:
        pmf="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt",
        est="diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.est.txt"
    benchmark:
        "bench/diffexp/{dataset}.{experiment1}.{group1}-vs-{experiment2}.{group2}.{settings}.txt"
    threads: 8
    shell:
        "{merfishtools} diffexp -t {threads} --pmf {output.pmf} {input} "
        "> {output.est}"


rule plot_expression_pmf:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.txt",
        raw_counts="counts/{dataset}.{experiment}.{group}.txt"
    output:
        "results/{context}/expression_pmf/{dataset}.{experiment}.{group}.{gene}.{settings}.expression_pmf.svg"
    script:
        "scripts/plot-expression-pmf.py"


rule plot_foldchange_cdf:
    input:
        pmf="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.txt",
        est="diffexp/{dataset}.{experiment}.{group1}-vs-{experiment}.{group2}.{settings}.est.txt"
    output:
        "results/{context}/foldchange_cdf/{dataset}.{experiment}.{group1}-vs-{group2}.{gene}.{settings}.foldchange_cdf.svg"
    script:
        "scripts/plot-foldchange-cdf.py"


rule expression_matrix:
    input:
        "expressions/{dataset}.{experiment}.{group}.{settings}.est.txt"
    output:
        "expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    script:
        "scripts/expression-matrix.py"


rule normalize_expression_matrix:
    input:
        expr="expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt",
        scales="normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    output:
        "normalized_expressions/{dataset}.{experiment}.{group}.{settings}.matrix.txt"
    run:
        expr = pd.read_table(input.expr, index_col=0)
        scales = pd.read_table(input.scales, index_col=0, squeeze=True, header=None)
        print(scales)
        expr *= scales[int(wildcards.experiment)]
        expr.to_csv(output[0], sep="\t")


rule count_matrix:
    input:
        "counts/{dataset}.{experiment}.{group}.txt"
    output:
        "counts/{dataset}.{experiment}.{group}.matrix.txt"
    script:
        "scripts/count-matrix.py"


rule plot_qq:
    input:
        lambda wildcards: matrices(wildcards.dataset, type=wildcards.type, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{type}.{settings}.qqplot.svg"
    params:
        experiments=lambda wildcards: experiments(wildcards.dataset)
    script:
        "scripts/plot-qq.py"


rule scale_factors:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "normalized_expressions/{dataset}.{settings}.scale_factors.txt"
    params:
        experiments=lambda wildcards: experiments(wildcards.dataset)
    script:
        "scripts/normalize.py"


rule plot_expression_dist:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.expression_dist.svg"
    script:
        "scripts/plot-expression-dist.py"


rule plot_overdispersion:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.overdispersion.svg"
    script:
        "scripts/plot-overdispersion.py"


rule plot_correlation:
    input:
        lambda wildcards: matrices(wildcards.dataset, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{settings}.correlation.svg"
    script:
        "scripts/plot-correlation.py"


rule plot_pca:
    input:
        lambda wildcards: matrices(wildcards.dataset, type=wildcards.type, settings=wildcards.settings)
    output:
        "results/{context}/{dataset}.{type}.{settings}.pca.svg"
    script:
        "scripts/plot-pca.py"


rule analyze_codebook:
    input:
        "data/{codebook}.codebook.txt"
    output:
        "neighbors/{codebook}.neighbors.txt"
    params:
        neighbor_dist=4#lambda wildcards: config["datasets"][wildcards.codebook.split(".")[0] + "Data"]["dist"]
    script:
        "scripts/codebook-neighbors.py"


rule plot_neighbor_bias:
    input:
        neighbors="neighbors/1001genes.1.neighbors.txt",
        counts=matrices("1001genesData", type="counts")
    output:
        "results/{context}/1001genesData.{settings}.neighbor_bias.svg"
    script:
        "scripts/plot-neighbor-bias.py"


rule simulate:
    input:
        mhd4="data/140genes.1.codebook.txt",
        mhd2="data/1001genes.1.codebook.txt"
    output:
        sim_counts_mhd4="data/simulated-MHD4.{mean}.all.txt",
        sim_counts_mhd2="data/simulated-MHD2.{mean}.all.txt",
        known_counts="data/simulated.{mean}.known.txt"
    params:
        cell_count=10
    script:
        "scripts/simulate-counts.py"


means = list(range(5, 40, 5))


rule plot_simulation:
    input:
        posterior_counts=expand("expressions/simulated-MHD{{dist}}.{mean}.all.{{settings}}.est.txt", mean=means),
        raw_counts=expand("counts/simulated-MHD{{dist}}.{mean}.all.txt", mean=means),
        known_counts=expand("data/simulated.{mean}.known.txt", mean=means)
    output:
        violin="results/{context}/simulation-MHD{dist}/MHD{dist}.error.{settings}.svg",
        scatter="results/{context}/simulation-MHD{dist}/MHD{dist}.scatter.{settings}.svg"
    params:
        means=means
    script:
        "scripts/plot-simulation.py"


rule plot_dataset_correlation:
    input:
        small=matrices("140genesData"),
        large=matrices("1001genesData"),
        small_counts=matrices("140genesData", type="counts"),
        large_counts=matrices("1001genesData", type="counts"),
        neighbors="neighbors/1001genes.1.neighbors.txt"
    output:
        "results/{context}/{settings}.dataset_correlation.svg"
    script:
        "scripts/plot-dataset-correlation.py"


rule figure1:
    input:
        b="results/paper/expression_pmf/140genesData.1.cell0.COL5A1.default.expression_pmf.svg",
        a="results/paper/expression_pmf/140genesData.1.cell0.CKAP5.default.expression_pmf.svg",
        c="results/paper/foldchange_cdf/140genesData.1.A-vs-B.COL5A1.default.foldchange_cdf.svg"
    output:
        "figures/fig1.svg"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("24.8cm", "6cm")
        a = sg.fromfile(input.a).getroot()
        b = sg.fromfile(input.b).getroot()
        c = sg.fromfile(input.c).getroot()
        b.moveto(294, 0)
        c.moveto(588, 0)

        la = sg.TextElement(0,10, "a", size=12, weight="bold")
        lb = sg.TextElement(294,10, "b", size=12, weight="bold")
        lc = sg.TextElement(588,10, "c", size=12, weight="bold")

        fig.append([a, b, c, la, lb, lc])
        fig.save(output[0])


rule figure2:
    input:
        a="results/paper/140genesData.default.expression_dist.svg",
        b="results/paper/default.dataset_correlation.svg",
        c="results/paper/simulation-MHD4/MHD4.error.default.svg",
        d="results/paper/simulation-MHD2/MHD2.error.default.svg"
    output:
        "figures/fig2.svg"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("16.2cm", "12cm")
        a = sg.fromfile(input.a).getroot()
        b = sg.fromfile(input.b).getroot()
        c = sg.fromfile(input.c).getroot()
        d = sg.fromfile(input.d).getroot()
        b.moveto(288, 0)
        c.moveto(0, 200)
        d.moveto(288, 200)

        la = sg.TextElement(0,10, "a", size=12, weight="bold")
        lb = sg.TextElement(288,10, "b", size=12, weight="bold")
        lc = sg.TextElement(0,210, "c", size=12, weight="bold")
        ld = sg.TextElement(288,210, "d", size=12, weight="bold")

        fig.append([a, b, c, d, la, lb, lc, ld])
        fig.save(output[0])


rule figure3:
    input:
        a="results/paper/140genesData.default.pca.svg",
        b="results/paper/140genesData.default.overdispersion.svg"
    output:
        "figures/fig3.svg"
    run:
        import svgutils.transform as sg
        fig = sg.SVGFigure("17.2cm", "11.5cm")
        a = sg.fromfile(input.a).getroot()
        b = sg.fromfile(input.b).getroot()
        b.moveto(0, 195)

        la = sg.TextElement(0,10, "a", size=12, weight="bold")
        lb = sg.TextElement(0,205, "b", size=12, weight="bold")

        fig.append([a, b, la, lb])
        fig.save(output[0])


rule svg2pdf:
    input:
        "{prefix}.svg"
    output:
        "{prefix}.pdf"
    shell:
        "rsvg-convert -f pdf {input} > {output}"
