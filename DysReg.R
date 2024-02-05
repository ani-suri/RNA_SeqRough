## Installing Packages 
library(tidyverse)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(oligo)
library(edgeR)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(ggpubr)
library(pd.hg.u219)
library(pd.hugene.2.0.st)
library(hugene20sttranscriptcluster.db)


#install.packages("BiocManager")
#install.packages("forcats")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("survminer")
#BiocManager::install("limma")
#BiocManager::install("pheatmap")
#BiocManager::install("org.Hs.eg.db")"

# Pulling the data using GSE ID 
# client specs 
GSE_ID = "GSE57071"
gse = getGEO(GSE_ID)[[1]]
length(gse) #Number of seq platforms, can

###
pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

### Normalization check 
#client spec 
exprs(gse) = log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

### Inspecting the clinical vars 
library(dplyr)
sampleInfo = pData(gse)
sampleInfo = select(sampleInfo, source_name_ch1,geo_accession)
sampleInfo = rename(sampleInfo,group = source_name_ch1, patient_acc=geo_accession)
sampleInfo 

### Renaming for grouping 
sampleInfo[1:3,1] = 'NonTreated'
sampleInfo[4:6, 1]= "Treated"
sampleInfo 

library(limma)
design = model.matrix(~0+sampleInfo$group)
colnames(design) = c("NonTreated","Treated")
design

#Logic can be changed:
#"It has been demonstrated that our power to detect differential expression can be improved if we filter lowly-expressed genes prior to performing the analysis. Quite how one defines a gene being expressed may vary from experiment to experiment, so a cut-off that will work for all datasets is not feasible. Here we consider that aroudn 50% of our genes will not be expressed, and use the median expression level as a cut-off."


# Metrics 
summary(exprs(gse))
### Gene selction using metric (can be changed as per client)
cutoff = median(exprs(gse)) #median cut off 
cutoff 
expressed_above_CO = exprs(gse) > cutoff #only genes above cutoff are expressed 
expressed_above_CO
keep_genes = rowSums(expressed_above_CO) > 2 #identify the genes expressed in more than 2 samples 
keep_genes
table(keep_genes) #genes removed and retained 
gse = gse[keep_genes,] #rewrite gse to contain only the expressed genes as per metric checks
gse

#Model fit 
## fit the model to the data using limFit 
# The result of which is to estimate the expression level in each of the groups that we specified.
fit = lmFit(exprs(gse), design)
head(fit$coefficients)

### Contrast 
contrasts = makeContrasts(NonTreated - Treated, levels=design)
fit2 = contrasts.fit(fit, contrasts)
fit2 = eBayes(fit2)
topTable(fit2) #list 
fData(gse)


# Final table 
final_table = fData(gse)
final_table 
final_table = select(final_table, ID, RANGE_START, RANGE_STOP, gene_assignment)
final_table$gene_assignment = ifelse(final_table$gene_assignment == "---", "NA", final_table$gene_assignment)

#Gene list 
gene_list = final_table$gene_assignment
gene_list_ten= data.frame(head(gene_list))
gene_list_ten

write.csv(gene_list_ten, "export10gene.csv", row.names = FALSE)


##heat map?? 
#test 

## Use to top 20 genes for illustration

#topN = 20
##
#ids_of_interest = mutate(final_table, Rank = 1:n()) %>% 
#  filter(Rank < topN) %>% c
#  pull(gene_assignment)

#gene_matrix <- exprs(gse)[ids_of_interest,]
