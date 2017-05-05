#!/usr/bin/evn Rscript
source('https://aidistan.github.io/gist/R/use.packages.R')
use.packages('parallel', 'igraph', 'ggplot2')

#
# Load inputs
#

cluster <- makePSOCKcluster(2)

# Phenotype-phenotyp similarities
phenotype_similarities <- read.table('../tmp/inner_phenotype_similarity.txt')

# Protein-protein interactions
ppi <- read.table('../tmp/inner_ppi.txt')
ppi_net <- graph_from_edgelist(as.matrix(ppi), directed = F)
gene_distances <- distances(ppi_net)

# Phenotype-gene relationships (strings as keys)
phenotype_gene_relationships <- read.table('../tmp/inner_phenotype_gene_relation.txt')

# Frequently-used numbers
gene_num <- max(ppi)
phenotype_num <- length(phenotype_similarities)

#
# Calculate gene-phenotype closeness
#

gene2phenotype_closeness <- matrix(data = 0, nrow = gene_num, ncol = phenotype_num)
for (phenotype in 1:phenotype_num) {
  phenotype_genes <- phenotype_gene_relationships[phenotype_gene_relationships[,1] == phenotype, 2]
  
  if (is.null(phenotype_genes)) {
    gene2phenotype_closeness[,phenotype] <- 0
  } else if (length(phenotype_genes) == 1) {
    gene2phenotype_closeness[,phenotype] <- exp(-gene_distances[,phenotype_genes]^2)
  } else {
    gene2phenotype_closeness[,phenotype] <- parRapply(cluster, exp(-gene_distances[,phenotype_genes]^2), sum)
  }
}

#
# Leave-one-out test
#

save(phenotype_similarities, gene_distances, gene2phenotype_closeness, phenotype_gene_relationships, file = '.Rdata')

# Test parallelly
system.time(
leave_one_out_results <- parRapply(cluster, phenotype_gene_relationships, function (row) {
  load('.Rdata')
  
  gene <- row[2]
  phenotype <- row[1]
  
  phenotype_genes <- phenotype_gene_relationships[
    phenotype_gene_relationships[,1] == phenotype & phenotype_gene_relationships[,2] != gene,
    2
  ]
  
  if (length(phenotype_genes) == 0) {
    gene2phenotype_closeness[,phenotype] <- 0
  } else if (length(phenotype_genes) == 1) {
    gene2phenotype_closeness[,phenotype] <- exp(-gene_distances[,phenotype_genes]^2)
  } else {
    gene2phenotype_closeness[,phenotype] <- apply(exp(-gene_distances[,phenotype_genes]^2), 1, sum)
  }

  gene_score <- cor(phenotype_similarities[,phenotype], t(gene2phenotype_closeness))
  gene_score[is.na(gene_score)] <- 0
  return(sum(quantile(gene_score, probs = seq(0, 1, 1e-04)) > gene_score[gene]) * 1e-04)
})
)

#
# Clean-up & Plot
#

stopCluster(cluster)

ggplot(data.frame(x = leave_one_out_results)) + theme_bw() +
  stat_ecdf(aes(x), geom = 'line') +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F)
