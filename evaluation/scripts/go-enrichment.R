library("GOstats")
library("biomaRt")
library("GO.db")
library("org.Hs.eg.db")
library("mutoss")

est <- read.table(snakemake@input[[1]], header = TRUE, stringsAsFactors = FALSE)
print(head(est$feat))

# translate gene names to entrez ids
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
ids <- getBM(attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = est$feat, mart = ensembl)
rownames(ids) <- ids$hgnc_symbol
est$entrez <- ids[est$feat, "entrezgene"]

# define foreground genes
significant <- est$diff_fdr <= 0.05
foreground <- est[significant, ]

print(head(foreground$entrez))

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

by <- BY(goterms$Pvalue, 0.05)
goterms$adjPvalue <- by[["adjPValues"]]

write.table(goterms, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE, sep = "\t")
