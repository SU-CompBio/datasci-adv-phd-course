setwd("~/Documents/Teaching/datasci-phd-course-2")
meta = read.csv('data/meta.csv', sep=',', header=T, row.names = 1)

fil = meta$Cell == 'A549'
meta = meta[fil,]

data = read.csv('results.csv', sep=',', header=T, row.names = 1)

BiocManager::install("piano")
install.packages('msigdbr')

library(piano)
library(msigdbr)

gene_set_kegg = msigdbr(species = 'Homo sapiens', category = 'C2', subcategory = 'KEGG')
gene_set_kegg = as.data.frame(gene_set_kegg[, c('gene_symbol', 'gs_name')])

gene_set_go = msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP')
gene_set_go = as.data.frame(gene_set_go[, c('gene_symbol', 'gs_name')])

fil = gene_set_kegg$gene_symbol %in% rownames(data)
gene_set_kegg = gene_set_kegg[fil, ]
fil = gene_set_go$gene_symbol %in% rownames(data)
gene_set_go = gene_set_go[fil, ]

gene_set_go = loadGSC(gene_set_go)
gene_set_kegg = loadGSC(gene_set_kegg)

gene_level_stat = data$log2FoldChange
names(gene_level_stat) = rownames(data)

kegg_results = runGSA(geneLevelStats = gene_level_stat, #ranking of genes 
                      geneSetStat = 'fgsea' , #method
                      gsc = gene_set_kegg,    #gene set
                      gsSizeLim = c(5, 300))
# https://www.pnas.org/content/102/43/15545

kegg_results =GSAsummaryTable(kegg_results)
head(kegg_results[order(kegg_results$`p adj (dist.dir.up)`),])
head(kegg_results[order(kegg_results$`p adj (dist.dir.dn)`),])

go_results = runGSA(geneLevelStats = gene_level_stat, #ranking of genes 
                      geneSetStat = 'fgsea' , #method
                      gsc = gene_set_go,    #gene set
                      gsSizeLim = c(5, 300))
go_results =GSAsummaryTable(go_results)

head(go_results[order(go_results$`p adj (dist.dir.up)`),])
head(go_results[order(go_results$`p adj (dist.dir.dn)`),])

BiocManager::install("dorothea")
BiocManager::install("viper")
library(dorothea)
library(viper)
library(dplyr)

data(dorothea_hs, package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

dset = data$log2FoldChange
names(dset) = rownames(data)
dset = as.data.frame(dset)

tf_activities <- run_viper(dset, regulons, 
                           options =  list(minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))
tf_activities[order(tf_activities),]


#### HW:
### Deseq2 A549 SARS-CoV-2 vs. Mock
### DEseq2 results A549 GO, KEGG, Dorothea + some other gene set enrichment
### Dorothea: scatter plot calu3 - a549 TF activity
