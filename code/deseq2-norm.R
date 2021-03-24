setwd("~/Documents/Teaching/datasci-phd-course-2")
library(DESeq2)

data = read.csv('data/GSE147507_RawReadCounts_Human.csv',
                sep=',', header=T, row.names = 1)
meta = read.csv('data/meta.csv',
                sep=',', header=T, row.names = 1)

dds = DESeqDataSetFromMatrix(countData = data,
                             colData = meta,
                             design = ~ Treatment)
vsd <- vst(dds, blind=TRUE)
write.csv(assay(vsd), 'norm_data.csv')
