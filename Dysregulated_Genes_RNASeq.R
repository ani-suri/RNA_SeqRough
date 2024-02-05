
  
  ```{r include=FALSE}
library(tidyverse)
library(magrittr)
library(GEOquery)
library(AnnotationDbi)
library(oligo)
library(edgeR)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(ggpubr)
```

```{r}
# Import file
gse <- getGEO("GSE193118")[[1]]
sampleInfo <- pData(gse)
''' table(sampleInfo$characteristics_ch1)
sampleInfo %<>% mutate(caco = ifelse(str_detect(characteristics_ch1,'healthy'),0,1))
sampleInfo %<>% mutate(id = rownames(sampleInfo))
sampleInfo %<>% mutate(age = `age:ch1`)
sampleInfo$age = as.numeric(sampleInfo$age)
hist(sampleInfo$age)
# age group <30, >30 
sampleInfo %<>% mutate(age30 = ifelse(age>=30, 1,0))
sampleInfo$caco = as.factor(sampleInfo$caco)
sampleInfo %>% ggplot(aes(age, colour=caco, fill=caco)) + geom_histogram(alpha=0.5)
sampleInfo %>% group_by(caco) %>% summarise(mean = mean(age), sd = sd(age))
sampleInfo %$% table(age30,caco)

# age group x>30, 30<x<60, 60<x
# limma does not covert it to dimmy
# need to create dummy
# sampleInfo %<>% mutate(age_young = ifelse(age < 30, 1, 0))
sampleInfo %<>% mutate(age_g = ifelse(is.na(age), NA,
                                      ifelse(age<=30, 0,
                                             ifelse(age>30 & age<=60, 1,
                                                    ifelse(age>=60, 2, NA)))))
sampleInfo %<>% mutate(age_mid = ifelse(age > 30 & age<=60, 1, 0))
sampleInfo %<>% mutate(age_old = ifelse(age >=60, 1, 0))'''
```


```{r}
# limma-voom
# raw counts
ex0 = read_csv("GSE193118_Readcount.csv")

# need to convert emsemble to entrez because GO analysis need entrezID
annotable = AnnotationDbi::select(org.Hs.eg.db, keys = ex0$ensembl_gene_id,
                                  columns = c("ENTREZID",'ENSEMBL','SYMBOL'),
                                  keytype = "ENSEMBL")
ex1 = merge(annotable[,c("ENTREZID",'ENSEMBL')], ex0, all.x = T, by.x = "ENSEMBL", by.y='ensembl_gene_id')

# entrez + samples
ex = ex1[,!names(ex1)%in% c('ENSEMBL')]
names(ex) = c('ENTREZID', sampleInfo$geo_accession)

# remove NA in entrez
ex3 = ex[!is.na(ex$ENTREZID),]
table(duplicated(ex3$ENTREZID))

# for duplicated entrez, take median for each gene
ex4 = ex3 %>% group_by(ENTREZID) %>% summarise_all(median)

# gene only and rowname = entrez
genes = ex4 %>% dplyr::select(sampleInfo$geo_accession)
rownames(genes) = ex4$ENTREZID

# quality control keep cpm>1 over 3 samples
cpm = cpm(genes)
keep = rowSums(cpm>1.5)>=3
kept_genes = genes[keep,]
gene_names_kept = ex4$ENTREZID[keep]
rownames(kept_genes) = gene_names_kept
dim(kept_genes)

# save rownames = entrez
saveRDS(kept_genes, 'kept_genes_185263.rds')
```


```{r}
# limma-voom
kept_genes = readRDS('kept_genes_185263.rds')
# https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/doc/usersguide.pdf 
# p46-48
# https://bioc.ism.ac.jp/packages/3.7/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

design = model.matrix(~0 + sampleInfo$caco + 
                        sampleInfo$caco:sampleInfo$age_mid + sampleInfo$caco:sampleInfo$age_old)

colnames(design) = c('con','case','con*mid','case*mid','con*old','case*old')
v = voom(kept_genes, design, plot = T)

vfit = lmFit(v, design)
# interaction assessment
vfit2 = contrasts.fit(vfit, c(0, 0, -1, 1, 0,  0))
# compare young vs mid
# age effect in sepsis   
# case mid     0 1 0 1 0 0 
# case young   0 1 0 0 0 0 
# effect       0 0 0 1 0 0 
# age effect in control
# con mid      1 0 1 0 0 0 
# con young    1 0 0 0 0 0 
# effect       0 0 1 0 0 0
# EMM additive effect
#              0 0 -1 1 0 0

# compare mid vs old
# age effect in sepsis   
# case old      0 1 0 0 0 1
# case mid      0 1 0 1 0 0 
# effect        0 0 0 -1 0 1
# age effect in control
# con old       1 0 0 0 1 0 
# con mid       1 0 1 0 0 0
# effect        0 0 -1 0 1 0
# EMM additive effect
#             0 0 1 -1 -1  1

# compare young vs old
# age effect in sepsis   
# case old        0 1 0 0 0 1
# case young      0 1 0 0 0 0  
# effect          0 0 0 0 0 1
# age effect in control
# con old       1 0 0 0 1 0 
# con young     1 0 0 0 0 0   
# effect        0 0 0 0 1 0 
# EMM additive effect
#             0 0 0 0 -1  1

efit = eBayes(vfit2)
top.table = topTable(efit, number = Inf)

top.table%<>% 
  mutate(ENTREZID = rownames(top.table)) %>%
  mutate(Significant = adj.P.Val < 0.05 &  abs(logFC) > 0.5 ) %>% 
  mutate(minuslog10=-log10(adj.P.Val)) %<>% 
  mutate(direction=ifelse(logFC>0, "up","down"))

# saveRDS(top.table,'top.table_185263.rds')
```

```{r}
# top.table = readRDS('top.table_185263.rds')
top.table %$% table(Significant)
top.table %>% filter(Significant) %$% table(direction)

# young vs mid
# Show in New Window
# Significant
# FALSE  TRUE 
# 12792  2680 
# direction
# down   up 
# 2044  636 

# mid vs old
# 0

# young vs old
# 0
```


```{r}
id = 10424
gexprs = kept_genes[rownames(kept_genes)==id, ]
df = data.frame(sampleInfo, geneA = t(gexprs))
df %>% group_by(caco) %>% ggplot(aes(x = as.factor(age_g), y=geneA, color = caco)) + geom_boxplot()
```




```{r}
# Volcano plot for soi
# create -log10(adj.Pval)
#top.table = readRDS('top.table_185263.rds')

v = top.table %>% 
  mutate(color_label=ifelse(Significant==TRUE & direction=="up","sig_up",
                            ifelse(Significant==TRUE & direction=="down","sig_down",
                                   ifelse(Significant==FALSE & 
                                            direction=="up","nonsig_up","nonsig_down")))) %>%
  ggplot(aes(x = logFC, y = minuslog10, col=color_label)) + 
  geom_point(size=1)+ 
  scale_color_manual(values = c("black","black", "#0066CC","#FF3300"))+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=14))+
  xlab("log2(FC)") +
  ylab("-log10(p.adjust)")+
  geom_hline(yintercept =-log10(0.05),linetype="dashed" ,color="grey")+
  geom_vline(xintercept = c(0.5,-0.5),linetype="dashed" ,color="grey")
#png('fig_volplot_185263.png', height = 7, width = 10, units='in', res= 400)
v
#dev.off()
```



```{r}
## GO enrichment analysis
top.table = readRDS('top.table_185263.rds')
sig_gene = top.table %>% filter(Significant)

direction<-c("up","down")
ontology<-c("BP","CC","MF")
ego_all<-data.frame()
for (i in direction){
  df<-sig_gene %>% filter(direction==i)
  for(j in ontology){
    ego <- enrichGO(gene = df[[11]], # entrez is 11th column
                    keyType = "ENTREZID",
                    OrgDb = org.Hs.eg.db, 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 1000,
                    ont = j,
                    readable = TRUE)
    ego<- data.frame(ego)
    ego <- ego %>% mutate(direction=i)
    ego <- ego %>% mutate(ontology=j)
    ego_all<-rbind(ego_all,ego)
  }}
saveRDS(ego_all,"ego_all_185263.rds")
```


```{r}
# go counts
ego_all<-readRDS("ego_all_185263.rds")
table(ego_all$direction,ego_all$ontology)
```

```{r}
# ggplot of GO
ego_all<-readRDS("ego_all_185263.rds")
ego_all %<>% mutate(minuslog_adjp=-log10(p.adjust)) 

# y axis will be automatically sorted by alphabatical order
ego_all$Description<- 
  factor(ego_all$Description,
         levels=ego_all[order(ego_all$minuslog_adjp,
                              ego_all$Description,decreasing = F),]$Description)
table(ego_all$direction, ego_all$ontology)
#        BP  CC  MF
#  down 500 122 125
#  up   252  35  12

# taka top 20 from down and up respectively
ego_all<-ego_all[order(ego_all$direction,ego_all$ontology,
                       ego_all$p.adjust,ego_all$Description),]
ego_all %<>% mutate(dir_new= 
                      ifelse(direction=="down","down-regulation", "up-regulation"))
up = ego_all %>% filter(direction=='up')
down = ego_all %>% filter(direction=='down')
ego_plot = rbind(up[1:20,], down[1:20,])

# ggplot
direction<-c("up","down")

for (i in direction){
  if (i=='up'){
    color_ = "#FF3300"
  } else{
    color_ = "#0066CC"  
  }
  gg = ego_plot %>% filter(direction ==i) %>% 
    ggplot(aes(x=minuslog_adjp,y=Description,
               shape=ontology,color=dir_new,size=Count)) +
    geom_point()+
    scale_color_manual(values = color_)+
    theme_light() +
    xlab("-Log10(p.adjust)") +
    ylab("")+
    labs(color="Direction", shape="Ontology",size="Gene count")+
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  name = assign(paste0('gg_',i),gg)
}
png('fig_go_plot.png', height=7, width=12, units= 'in', res =400)
ggarrange(gg_down,gg_up)
dev.off()
```


## KEGG pathway enrichment analysis
```{r}
top.table = readRDS('top.table_185263.rds')
sig_gene = top.table %>% filter(Significant)

KEGG_all<-data.frame()
for (i in direction){
  df = sig_gene %>% filter(direction==i)
  k = enrichKEGG(gene = df[[11]], organism = 'hsa',pvalueCutoff = 0.05,)
  k = as.data.frame(k)
  k %<>% mutate(direction=i)
  KEGG_all = rbind(KEGG_all,k)
}
saveRDS(KEGG_all,'kegg_all_185263.rds')
```


```{r}
# linear model
n= 58
m = ncol(df)
result = data.frame()
for(i in seq(n+1,m)){
  gene = names(df)[i]
  fit = lm(df[[i]] ~ caco + age + caco*age, data=df)
  coef = summary(fit)$coefficients
  result[i-n,1] = gene
  result[i-n,2] = coef[2,1]
  result[i-n,3] = coef[3,1]
  result[i-n,4] = coef[4,1]
  result[i-n,5] = coef[2,4]
  result[i-n,6] = coef[3,4]
  result[i-n,7] = coef[4,4]
}
names(result) = c('gene','coef_caco','coef_age','coef_caco_age',
                  'p_valuecaco','p_valueage','p_valueinter')
result %<>% mutate(fdr = p.adjust(p_valueinter, method = 'fdr',n=m-n))
sig = result %>% filter(fdr < 0.05)
sig
```







