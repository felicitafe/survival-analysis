library(TCGAbiolinks)
library(TCGAbiolinksGUI.data)
library(limma)
library(edgeR)
#library(glmnet)
library(factoextra)
library(FactoMineR)
library(caret)
library(SummarizedExperiment)
library(gplots)
library(survival)
library(survminer)
library(RColorBrewer)
#library(gProfileR)
library(genefilter)
library(BiocManager)
library(data.table)
library(tidyverse)
library(GenomicFeatures)
library(biomaRt)
library(DESeq2)
library(ggfortify)
library(ggrepel)
library(Glimma)
library(apeglm)
library(ashr)
library(circlize)
library(openxlsx)
library(GSVA)
library(tximport)
library(ggplot2)

GDCprojects = getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

# Raw counts  with HTSeq estimation
query_TCGA_BRCA = GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling", 
                            workflow.type = "STAR - Counts", data.type = "Gene Expression Quantification")


BRCA_res$sample_type <- as.factor(BRCA_res$sample_type)
summary(BRCA_res$sample_type)

# Download the files from the query, creates a new folder
GDCdownload(query = query_TCGA_BRCA)
tcga_BRCA_data = GDCprepare(query_TCGA_BRCA)


clinical = tcga_BRCA_data@colData
countdata = assay(tcga_BRCA_data)

# retrieve only the relevant information for now
sampledata <- data.frame(clinical$barcode, clinical$days_to_death, clinical$patient , clinical$sample, clinical$sample_type_id , clinical$sample_id, clinical$sample_type, clinical$primary_diagnosis, clinical$vital_status, clinical$paper_BRCA_Subtype_PAM50, clinical$paper_vital_status, clinical$paper_days_to_birth, clinical$paper_days_to_death, clinical$paper_days_to_last_followup, clinical$paper_age_at_initial_pathologic_diagnosis, clinical$paper_pathologic_stage, clinical$days_to_diagnosis, clinical$last_known_disease_status, clinical$tissue_or_organ_of_origin, clinical$days_to_last_follow_up, clinical$age_at_diagnosis, clinical$primary_diagnosis)

clin_df = clinical[clinical$definition == "Primary solid Tumor",c("patient",
                     "vital_status", #whether the patient is alive or dead
                     "days_to_death", #the number of days passed from initial diagnosis to the death
                     "days_to_last_follow_up", #the number of days passed from initial diagnosis to last visit
                     "paper_BRCA_Subtype_PAM50")]




keep <- rowSums(countdata) > 10
countdata <- countdata[keep,]
sampledata$clinical.sample_type <- as.factor(sampledata$clinical.sample_type)
sampledata$clinical.vital_status <- as.factor(sampledata$clinical.vital_status)
sampledata$clinical.sample_type<-relevel(sampledata$clinical.sample_type,ref="Solid Tissue Normal")




sampledatasubtype <- sampledata[sampledata$clinical.paper_BRCA_Subtype_PAM50 %in% c("LumA", "Normal"), ]
sampledatalum <- sampledata[sampledata$clinical.paper_BRCA_Subtype_PAM50 %in% c("LumA"), ]
clin_df_sub <- clin_df[clin_df$paper_BRCA_Subtype_PAM50 %in% c("LumA", "Normal"), ]


countdatasub <- countdata[, sampledatasubtype$clinical.barcode]

# generate a matrix with counts for the barcodes in sampledatasubtype, and group according to the subtype
dds <- DESeqDataSetFromMatrix(countdatasub, sampledatasubtype, ~clinical.paper_BRCA_Subtype_PAM50)


ddsObjI <- DESeq(dds)
normcounts <- counts(ddsObjI, normalized=T)
vstcounts <- vst(ddsObjI, blind=TRUE)
counts <- vstcounts@assays@data@listData[[1]]


bxgenes <- c("ENSG00000143409.15", "ENSG00000133069.17")
hub <- c("FAM63A", "TMCC2" )

exp.mat <- counts[paste(bxgenes), ]
print(exp.mat)

rownames(exp.mat) <- c("FAM63A", "TMCC2")



res1 <- results(ddsObjI, alpha=0.05)
resultsNames(ddsObjI)
results1 <- as.data.frame(res1)
results1 <- res1[paste(bxgenes), ]

results1 <- as.data.frame(results1) %>% rownames_to_column("gene_id")
#results1$padj <- round(results1$padj,3)


clin_df_sub$deceased = clin_df_sub$vital_status == "Dead"

clin_df_sub$overall_survival = ifelse(clin_df_sub$deceased,
                                  clin_df_sub$days_to_death,
                                  clin_df_sub$days_to_last_follow_up)

# Surv(clin_df$overall_survival, clin_df$deceased)

fit = survfit(Surv(overall_survival, deceased) ~paper_BRCA_Subtype_PAM50, data=clin_df_sub)
print(fit)

pval = surv_pvalue(fit, data=clin_df_sub)$pval



ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
ourFilterType <- "ensembl_gene_id_version"
attributeNames <- c('ensembl_gene_id_version','entrezgene_accession')


# run the query
filterValues <- rownames(normcounts)
annot <- getBM(attributes=attributeNames, 
               filters = ourFilterType, 
               values = filterValues, 
               mart = ensembl)

annot <- annot[!duplicated(annot$entrezgene_accession),]

normcounts2 <- as.data.frame(normcounts, row.names=NULL)%>%
  rownames_to_column("ensembl_gene_id_version") %>% 
  left_join(annot, c("ensembl_gene_id_version"))%>%
  drop_na("entrezgene_accession") %>%
 column_to_rownames("entrezgene_accession")


# ensembl id leri sil, sadece gen isimleri ve hasta datasını sakla
normcounts3 <- normcounts2[,-1]


# divide patients in two groups, up and down regulated. If the patient expression is greater or equal to the median we put it
# among the "up-regulated", otherwise "down-regulated"
clin_df_sub$gene = as.factor(ifelse(clin_df_sub$gene_value >= median, "UP", "DOWN"))
# fit a gene survival model
fit = survfit(Surv(overall_survival, deceased) ~gene, data=clin_df_sub)
# extract the survival p-value
pval = surv_pvalue(fit, data=clin_df_sub)$pval
# Kaplan-Meier plot
pdf("survplot_plots")
plotsurv <- ggsurvplot(fit, data=clin_df_sub, pval=T, risk.table=T, title=paste(rownames(exp.mat)[1]))
plot_list[[i]] = plotsurv
dev.off()
# clin_df_sub <- clin_df_sub[,-c(25,26)]
}
