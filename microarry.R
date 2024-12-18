getwd()
library(GEOquery)
library(limma)
BiocManager::install("umap")
library(umap)
library("oligo")
BiocManager::install("Biobase")
library("Biobase")
library(ggplot2)

raw_data_dir <- "D:/data_need/Microarry_data"
SDRF <-  read.csv("SDRF.csv", sep = ",")
rownames(SDRF) <- SDRF$Arrat.Data.file
SDRF <- AnnotatedDataFrame(SDRF)

#Reading CEL File

raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, SDRF$Arrat.Data.file),
                                 verbose = FALSE, phenoData = SDRF)

stopifnot(validObject(raw_data))

head(Biobase::pData(raw_data))
names(Biobase::pData(raw_data))

##Quality control of raw data

Biobase::exprs(raw_data)[1:5,1:5]

exp_raw <- log2(Biobase::exprs(raw_data))

# Box Plot of the raw data exoression data

BOXPLOT_Raw <- oligo::boxplot(raw_data, target = "core",las=2,
               main = " Boxplot of log2-intensitites for the raw data")

# constructing pca plot of rawdataexpression

PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentvar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentvar[2]/percentvar[1])
dataGG <- data.frame(PC1 =PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     phenotype = pData(raw_data)$expression.type)

ggplot(dataGG, aes(PC1,PC2)) +
  geom_point(aes(shape = phenotype, colour = phenotype)) +
  ggtitle("PCA PLOT OF THE LOG-TRANSFORMED RAW EXPRESSION DATA") +
  xlab(paste0("PC1, Varexp: ", percentvar[1],"%"))+
  ylab(paste0("PC2, Varexp: ", percentvar[2],"%"))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio)+
  scale_shape_manual(values= c(4, 15))+
  scale_colour_manual(values= c("red", "blue"))

## rma calibration of the data (robust multi-array average) (quantile normalization)

palmieri_eset_norm <- oligo::rma(raw_data)


# normalize data constructing pca plot of rawdataexpression

exp_palmieri <- Biobase::exprs(palmieri_eset_norm)

PCA_raw <- prcomp(t(exp_palmieri), scale. = FALSE)

percentvar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentvar[2]/percentvar[1])
dataGG <- data.frame(PC1 =PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     phenotype = pData(palmieri_eset_norm)$expression.type)

ggplot(dataGG, aes(PC1,PC2)) +
  geom_point(aes(shape = phenotype, colour = phenotype)) +
  ggtitle("PCA PLOT OF THE LOG-TRANSFORMED RAW EXPRESSION DATA") +
  xlab(paste0("PC1, Varexp: ", percentvar[1],"%"))+
  ylab(paste0("PC2, Varexp: ", percentvar[2],"%"))+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio)+
  scale_shape_manual(values= c(4, 15))+
  scale_colour_manual(values= c("red", "blue"))

## box plot for normalize data

pdf("")
BOXPLOT_Normalized <- oligo::boxplot(palmieri_eset_norm, target = "core",
                                     las=2,
                              main = " Boxplot of log2-intensitites for the normalized data")

## heat map construction ##
row.names(pData(palmieri_eset_norm))
library(stringr)
library(pheatmap)
library(RcolorBrewer)
install.packages("RcolorBrewer")
BiocManager::install("RcolorBrewer")

pData(palmieri_eset_norm)

Cell_type <- ifelse(str_detect(pData(palmieri_eset_norm)$expression.type ,
                               "Exhausted CD8"), "Exhausted CD8","Effector CD8")

Annotation_for_heatmap <- data.frame(celltype = Cell_type)

Annotation_for_heatmap

row.names(Annotation_for_heatmap) <- row.names(pData(palmieri_eset_norm))
row.names(pData(palmieri_eset_norm))
row.names(Annotation_for_heatmap)

dists <- as.matrix(dist(t(exp_palmieri), method = "manhattan"))

rownames(dists) <- row.names(pData(palmieri_eset_norm))


hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255))

colnames(dists) <- NULL
diag(dists) <- NA

ann_colurs <- list(
  cell_type = c("Exhausted CD8" = "green", "Effector CD8" = "red" )
)
pheatmap(dists, col = (hmcol),
         annotation_row = Annotation_for_heatmap,
         annotation_colors = ann_colurs,
         legend = TRUE,
         treeheight_row = 0,
         legend_breaks = c(min(dists, na.rm =TRUE),
                           max(dists, na.rm =TRUE)),
         legend_labels = (c("small distance", "large distance")),
         msin = "clustering heatmap")


## probe annotations ###

BiocManager::install("mouse4302.db")
library(mouse4302.db)
library(dplyr)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)



###  annotation with mouse4302.db ##

Annotations <- AnnotationDbi::select(mouse4302.db,
                                     keys = (featureNames(palmieri_eset_norm)),
                                     columns = c("SYMBOL", "GENENAME"),
                                     keytype = "PROBEID")

Annotations <- subset(Annotations, !is.na(SYMBOL))

## removing multiple mapping ###


anno_groped <- group_by(Annotations, PROBEID)
anno_summarized <- 
  dplyr::summarize(anno_groped, no_of_matches = n_distinct(SYMBOL))

anno_filtered <- filter(anno_summarized, no_of_matches > 1)

probe_stats <- anno_filtered

nrow(probe_stats)

ids_to_exlude <- (featureNames(palmieri_eset_norm) %in% probe_stats$PROBEID)

table(ids_to_exlude)

Final_counts <- subset(palmieri_eset_norm, !ids_to_exlude)

validObject(Final_counts)

fData(Final_counts)$PROBEID <- rownames(fData(Final_counts))

fData(Final_counts) <- left_join(fData(Final_counts), Annotations, by = "PROBEID")

### restore rownames after left_join###

rownames(fData(Final_counts)) <- fData(Final_counts)$PROBEID
validObject(Final_counts)

Biobase::pData(Final_counts)
fData(Final_counts) <- fData(Final_counts)[, !names(fData(Final_counts)) %in% c("SYMBOL.y", "GENENAME.y")]

fData(Final_counts) <- fData(Final_counts) %>%
  dplyr::rename(SYMBOL = SYMBOL.x, GENENAME = GENENAME.x)


write.csv(Final_counts, "Final_counts.csv")

write.csv(fData(Final_counts)[, c("PROBEID", "SYMBOL", "GENENAME")], "Final_counts_1.csv")

####till linear model ####



individual <-
  as.character(Biobase::pData(Final_counts)$expression.type)


Cell_type_2 <- str_replace_all(Biobase::pData(Final_counts)$expression.type,
                               " ", " ")

Cell_type_2 <- ifelse(str_detect(pData(Final_counts)$expression.type ,
                                 "Effector CD8"), "Effector_CD8","Exhausted_CD8")

disign_cell_type <- model.matrix(~ 0 + Cell_type_2)

colnames(disign_cell_type)[1:2] <- c("Effector_CD8", "Exhausted_CD8")

rownames(disign_cell_type) <- individual

write.csv(disign_cell_type, "design matrix.csv")

### contrasts and hypothesis tests ##

Fit <- lmFit(disign_cell_type)

### white taking constrast, you should place first those samples ##

contrast_matrix <- makeContrasts("Effector_CD8-Exhausted_CD8", levels = disign_cell_type )

### Building Linear Models contrast design ##

celltype_fit <- eBayes(contrasts.fit(lmFit(Final_counts, disign_cell_type),contrast_matrix))

##Extracting final results ###

table <- topTable(celltype_fit, number = Inf)

write.csv(table, "top_table_final.csv")


### excluding genes with no symbol

table <- subset(table, !is.na(SYMBOL))

colnames(table)

## Most significant Differentially Expressed Genes

DEG_Genes <- subset(table, adj.P.Val < 0.05)

write.csv(DEG_Genes,"DEG_Genes_Mostsignificat.csv")

## Valcano Plot ##
EnhancedVolcano(table,
               lab = table$SYMBOL,
               x = "logFC",
               y = "P.Value",
               ylim = c(0, -log10(10e-12)),
               pCutoff = 0.05,
               FCcutoff = 0.5,
               title = "Effector_CD8 vs Exhausted_CD8")
str(table$p.value)














