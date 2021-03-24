setwd("~/Documents/Teaching/datasci-phd-course-2")
library(DESeq2)

data = read.csv('data/GSE147507_RawReadCounts_Human.csv',
                sep=',', header=T, row.names = 1)
meta = read.csv('data/meta.csv',
                sep=',', header=T, row.names = 1)

fil = meta$Cell == 'Calu3'
meta = meta[fil,]

data = data[, rownames(meta)]

gene_exp = apply(data, 1, sum)
fil = gene_exp > 0
data = data[fil,]

meta$Treatment = as.factor(meta$Treatment)
meta$Treatment = relevel(meta$Treatment, 'Mock')

dds = DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, 'results.csv')
