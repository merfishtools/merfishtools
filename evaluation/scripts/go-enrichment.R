library("GOstats")
library("biomaRt")
library("GO.db")
library("org.Hs.eg.db")
library("mutoss")

est <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)

# translate gene names to entrez ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = est$feat, mart = ensembl)
rownames(ids) <- ids$hgnc_symbol
est$entrez <- ids[est$feat, "entrezgene"]

# define foreground genes
significant <- est$diff_fdr <= 0.05
foreground <- est[significant, ]

# define test parameters
params <- new("GOHyperGParams",
              geneIds = foreground$entrez,
              universeGeneIds = est$entrez,
              ontology = "BP",
              annotation = "org.Hs.eg",
              pvalueCutoff = 0.05,
              conditional = TRUE,
              testDirection = "over")
results <- hyperGTest(params)
goterms <- summary(results)

# correct for multiple testing. We use Benjamini-Yekuteli here, because the performed tests are strongly dependent
by <- BY(goterms$Pvalue, 0.05)
goterms$adjPvalue <- by[["adjPValues"]]

# Compute the DAG of significant go terms
graph <- inducedTermGraph(results, id = goterms[goterms$adjPvalue <= 0.05, ]$GOBPID)

# Plot the DAG.
pdf(snakemake@output[["graph"]])
plotGOTermGraph(graph, results)
dev.off()


write.table(goterms, file = snakemake@output[["table"]], row.names = FALSE, quote = FALSE, sep = "\t")
