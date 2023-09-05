library(tximport)
library(dplyr)

files_dir <- "~/NeuroGenom/Final Project/Data"
# List of all the abundance files in the directory
data_files <- list.files(files_dir, pattern = "abundance_.*\\.tsv$", full.names = TRUE)
tx2gene=read.csv("~/NeuroGenom/mus_musculus/mus_musculus/transcripts_to_genes.txt",sep='\t',header = FALSE)
res1<-tximport(data_files,type="kallisto",tx2gene=tx2gene,ignoreAfterBar=TRUE)
counts<-data.frame(round(res1$counts))
names(counts)<-c("C1","C2","C3","KO1","KO2","KO3")
rownames(counts) <- rownames(res1$counts)
write.table(counts, file = "~/NeuroGenom/Final Project/Data/counts_df.txt", sep = "\t", row.names = TRUE, col.names = NA)


library("DESeq2")
cts<-counts # the counts file
coldata<-read.table(file="~/NeuroGenom/Final Project/conditions_file.txt", sep = '\t') # the conditions file
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Assign the gene symbols as row names in the counts table
# deal with "AC148174.4,AC132307.6" 
matched_gene_symbols <- tx2gene$V3[match(row.names(res)[1:(nrow(res)-1)], tx2gene$V2)]
matched_gene_symbols<-c(matched_gene_symbols,"Null")
res$gene_symbol<-matched_gene_symbols
write.table(res, file = "~/NeuroGenom/Final Project/Data/Deseq_res.txt", sep = "\t", row.names = TRUE, col.names = NA)


#analysing the data
alpha<-0.05
low_p<-subset(res,padj<alpha)

groupA<-subset(low_p,log2FoldChange>0)
groupB<-subset(low_p,log2FoldChange<0)
#for usind David database
write.table(groupA$gene_symbol, file = "~/NeuroGenom/Final Project/Data/groupA_gene_symbols.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(groupB$gene_symbol, file = "~/NeuroGenom/Final Project/Data/groupB_gene_symbols.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

sortedA<-groupA[order(-groupA$log2FoldChange),]
sortedB<-groupB[order(groupB$log2FoldChange),]



