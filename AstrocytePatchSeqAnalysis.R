#Script for processing Patch-Seq data from Allen et. al. 2021

library(Seurat)
library(ggplot2)
library(pheatmap)
library(viridis)

#Create seurat object from each library--------------
#Library 1---------
raw <- Read10X(data.dir = "/Solo.out/Gene/raw/")
Lib1_deep <- CreateSeuratObject(counts = raw, min.cells = 1)
Lib1_deep
Lib1_deep@meta.data

#Cell number
Lib1_deep <- AddMetaData(Lib1_deep, metadata = colnames(Lib1_deep), col.name = "CellNumber")
#Celltype
Lib1_deep <- AddMetaData (Lib1_deep, metadata = "DenseBulbousAstrocyte", col.name = "CellType")
Lib1_deep@meta.data$CellType[19] = "TotalRNA"
#Library
Lib1_deep@meta.data$Library = 1
#Individual
Individual <- c("GW18_210617", "GW18_210617", "GW18_210617", "GW18_210617",
                "GW24_210630", "GW24_210630", "GW24_210630",
                "GW21_210726", "GW21_210726", "GW21_210726", "GW21_210726",
                "GW20_210730", "GW20_210730",
                "GW19_210806", "GW19_210806", "GW19_210806", "GW19_210806", "GW19_210806",
                "TotalRNA")
names(Individual) <- Lib1_deep@meta.data$CellNumber
Lib1_deep@meta.data$Individual <- Individual
#Age
Age <- c("GW18", "GW18", "GW18", "GW18",
         "GW24", "GW24", "GW24",
         "GW21", "GW21", "GW21", "GW21",
         "GW20", "GW20",
         "GW19", "GW19", "GW19", "GW19", "GW19",
         "TotalRNA")
names(Age) <- Lib1_deep@meta.data$CellNumber
Lib1_deep@meta.data$Age <- Age

Lib1_deep@meta.data
#saveRDS(Lib1_deep, file = "Library1.rds")
#Lib1 <- readRDS(file = "Library1.rds")

#Library 2--------------------------
raw <- Read10X(data.dir = "/Solo.out/Gene/raw/")
Lib2_deep <- CreateSeuratObject(counts = raw, min.cells = 1)
Lib2_deep
Lib2_deep@meta.data

#Cell number
Lib2_deep <- AddMetaData(Lib2_deep, metadata = colnames(Lib2_deep), col.name = "CellNumber")
Lib2_deep@meta.data$CellNumber[4] <- "Unknown"
#CellType
CellType <- c("DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "Unknown", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte",
              "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "Neuron", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "Neuron", "Unknown",
              "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "Neuron", "Neuron", "DenseBulbousAstrocyte", "Neuron", "DenseBulbousAstrocyte", "DenseBulbousAstrocyte", "Neuron")
names(CellType) <- Lib2_deep@meta.data$CellNumber
Lib2_deep <- AddMetaData (Lib2_deep, metadata = CellType, col.name = "CellType")
Lib2_deep@meta.data$CellType[4] <- "Unknown"
Lib2_deep@meta.data$CellType <- factor(Lib2_deep@meta.data$CellTyp, levels = c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte", "Neuron", "Unknown"))
#Library
Lib2_deep@meta.data$Library <- 2
#Individual
Individual <- c("GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930", "GW19_210930",
                "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112", "GW24_211112",
                "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115", "GW19_211115")
names(Individual) <- Lib1_deep@meta.data$CellNumber
Lib2_deep@meta.data$Individual <- Individual
#Age
Age <- c("GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19",
         "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24", "GW24",
         "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19", "GW19")
names(Age) <- Lib2_deep@meta.data$CellNumber
Lib2_deep@meta.data$Age <- Age

Lib2_deep@meta.data
#saveRDS(Lib2_deep, file = "Library2.rds")
#Lib2 <- readRDS(file = "Library2.rds")

#Library3------
raw <- Read10X(data.dir = "/Solo.out/Gene/raw/")
Lib3_deep <- CreateSeuratObject(counts = raw, min.cells = 1)
Lib3_deep
Lib3_deep@meta.data

#Cell number
Lib3_deep <- AddMetaData(Lib3_deep, metadata = colnames(Lib3_deep), col.name = "CellNumber")
Lib3_deep@meta.data$CellNumber[29] = "Uknown"

#CellType
CellType <- c("DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "Neuron", "Neuron", "Neuron", "Neuron", "DenseSmoothAstrocyte", "Neuron", "Neuron", "Neuron",
              "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "Neuron", "DenseSmoothAstrocyte", "Neuron",
              "Neuron", "Neuron", "Neuron", "Neuron", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "DenseSmoothAstrocyte", "Neuron", "Neuron", 
              "Neuron", "Astrocyte", "Uknown")
names(CellType) <- Lib3_deep@meta.data$CellNumber
Lib3_deep <- AddMetaData (Lib3_deep, metadata = CellType, col.name = "CellType")
Lib3_deep@meta.data$CellType[29] = "Unknown"
Lib3_deep@meta.data$CellType <- factor(Lib3_deep@meta.data$CellTyp, levels = c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte", "Neuron", "Astrocyte", "Unknown"))
#Library
Lib3_deep@meta.data$Library <- 3
#Individual
Individual <- c("GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210", "GW21_211210",
                "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", 
                "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211", "GW23_211211")
names(Individual) <- Lib3_deep@meta.data$CellNumber
Lib3_deep@meta.data$Individual <- Individual
#Age
Age <- c("GW21", "GW21", "GW21", "GW21", "GW21", "GW21", "GW21", "GW21", "GW21", "GW21",
         "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", 
         "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23", "GW23")
names(Age) <- Lib3_deep@meta.data$CellNumber
Lib3_deep@meta.data$Age <- Age

Lib3_deep@meta.data
#saveRDS(Lib3_deep, file = "Library3.rds")
#Lib3 <- readRDS(file = "Library3.rds")

rm(raw)


#Merge the libraries together--------------
PS_deep <- merge(x=Lib1_deep, y=list(Lib2_deep, Lib3_deep))
PS_deep
PS_deep@meta.data
#saveRDS(PS_deep, "PS_deep.rds")
#PS_deep <- readRDS("PS_deep.rds")

#Process the data--------------
#Select just the cells
Idents(PS_deep) <- "CellType"
VlnPlot(PS_deep, features = "nCount_RNA", group.by = "CellType")
PS_cells_deep <- subset(PS_deep, idents = c("TotalRNA", "Unknown"), invert = T)

#Remove cells with incorrect morphology and degraded cDNA
#Bad morphology: C1, C2, C3, C4, C81
#Bad cDNA: C14, C15, C28, C33, C48
Idents(PS_cells_deep) <- "CellNumber"
PS_cells_deep <- subset(PS_cells_deep, idents = c("C1", "C2", "C3", "C4", "C14", "C15", "C28", "C33", "C48", "C81"), invert = T)

#Remove cells with <1000 reads per cell
#Idents(PS_cells_deep) <- "CellNumber"
#PS_high_deep <- subset(PS_cells_deep, idents = c("C5", "C6", "C8", "C9", "C12", "C13", "C17", "C18", "C21", "C24", "C46"), invert = T)

#Remove cells with <100,00 reads per cell
Idents(PS_cells_deep) <- "CellNumber"
PS_high_deep <- subset(PS_cells_deep, idents = c("C5", "C6", "C8", "C9", "C10", "C11", "C12", "C13", "C17", "C18", "C21", "C24", "C37", "C46"), invert = T)


#Process these cells
PS_high_deep <- NormalizeData(object = PS_high_deep)
PS_high_deep <- FindVariableFeatures(object = PS_high_deep)
PS_high_deep <- ScaleData(object = PS_high_deep, features = rownames(PS_high_deep))
PS_high_deep <- RunPCA(object = PS_high_deep, npcs = 30)
ElbowPlot(PS_high_deep)
PS_high_deep <- FindNeighbors(object = PS_high_deep)
PS_high_deep <- FindClusters(object = PS_high_deep)
PS_high_deep <- RunTSNE(object = PS_high_deep, dims = 1:7, perplexity=10)
PS_high_deep <- RunUMAP(object = PS_high_deep, dims = 1:7)
DimPlot(object = PS_high_deep, reduction = "tsne")
DimPlot(object = PS_high_deep, reduction = "umap")

#saveRDS(PS_high_deep, file = "PS_high_deep.rds")
#PS_high_deep <- readRDS(file = "PS_high_deep.rds")

#Refine cell types and remove cells with poor gene expression-------------

#Construct list of cell type marker genes
Astro_genes <- c("SOX9", "ALDH1L1", "GJA1", "AQP4", "GFAP")
RG_genes <- c("SOX2", "VIM", "FABP7", "CRYAB", "HOPX")
Oligo_genes <- c("OLIG2",  "PDGFRA", "SOX10", "MBP")
Glia_genes <- c(Astro_genes, RG_genes, Oligo_genes)

EN_genes <- c("EOMES", "NEUROD2", "NEUROD6", "SATB2")
IN_genes <- c("DLX1", "DLX2",   "DLX5")
Div_genes <- c("MKI67", "CENPF", "TOP2A")
MG_genes <- c("AIF1", "CSF3R", "CX3CR1")
Endo_genes <- c("FN1", "COL1A1")
All_genes <- c(Astro_genes, RG_genes, Oligo_genes, Div_genes, EN_genes, IN_genes, MG_genes, Endo_genes)

DoHeatmap(PS_high_deep, features = c(All_genes, "ANGPTL4", "ITGB4"), slot = "data", draw.lines = T, group.by = "CellType") + scale_fill_gradient(low="grey91", high="blue")

#Plot by cell type then cell number to identify cells with improper gene expression
Idents(PS_high_deep) <- "CellType"
DB <- subset(PS_high_deep, idents = "DenseBulbousAstrocyte")
DS <- subset(PS_high_deep, idents = "DenseSmoothAstrocyte")
Neu <- subset(PS_high_deep, idents = "Neuron")

DoHeatmap(DB, features = All_genes, slot = "data", draw.lines = F, group.by = "CellNumber") + scale_fill_gradient(low="grey91", high="blue")
DoHeatmap(DS, features = All_genes, slot = "data", draw.lines = F, group.by = "CellNumber") + scale_fill_gradient(low="grey91", high="blue")
DoHeatmap(Neu, features = All_genes, slot = "data", draw.lines = F, group.by = "CellNumber") + scale_fill_gradient(low="grey91", high="blue")

#Based on gene expression, remove cells with gene expression that does not match their cell type annotation
#C63 and C66 express microglial genes.
#C47, 50, 53 are neurons that express high levels of glial genes
#C54, C60, C19, C40, C44, C43, C78 express high levels of genes from multiple cell types
Idents(PS_high_deep) <- "CellNumber"
PS_high_deep <- subset(PS_high_deep, idents = c("C47", "C50", "C53", "C63", "C66", "C54", "C60", "C19", "C40", "C44", "C43", "C78"), invert = T)

#Also based on expression on heatmap, refine cell type names
PS_high_deep$CellType[c("C56", "C57", "C58", "C61", "C62", "C67", "C71", "C72")] <- "IN" #Interneurons
PS_high_deep$CellType[c("C32", "C59", "C69", "C70", "C73", "C79", "C80")] <- "EN" #Excitatory neurons
PS_high_deep$CellType[c("C30", "C75", "C76")] <- "DenseSmooth-Dividing" #Dividing

PS_high_deep$CellType <- factor(PS_high_deep$CellType, levels = c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte", "DenseSmooth-Dividing", "IN", "EN"))
DoHeatmap(PS_high_deep, features = All_genes, slot = "data", draw.lines = T, group.by = "CellType") + scale_fill_gradient(low="grey91", high="blue")

#saveRDS(PS_high_deep, file = "PS_high_deep.rds")
#PS_high_deep <- readRDS(file = "PS_high_deep.rds")

#Perform differential gene expression analysis-----------------
Idents(PS_high_deep) <- "CellType"
Markers <- FindAllMarkers(PS_high_deep, only.pos = T, logfc.threshold = 0.1)
View(Markers)
#write.table(Markers, sep = "\t", file = "PS_high_deep_ClusterMarkers.txt")
#Markers <- read.table(sep = "\t", file = "PS_deep_ClusterMarkers.txt")

#Plot on heatmap
DoHeatmap(PS_high_deep, features = c(Markers$gene), slot = "data", draw.lines = T, group.by = "CellType", raster = F) + scale_fill_gradient(low="grey91", high="blue")
plt <- subset(PS_high_deep, idents = c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte"))
DoHeatmap(plt, features = Markers$gene[Markers$cluster %in% c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte")], slot = "data", draw.lines = T, group.by = "CellType", raster = F) + scale_fill_gradient(low="grey91", high="blue")

#Filter to most specific astrocyte genes---------
Markers_filt <- Markers[Markers$cluster %in% c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte"),]

#Remove mitochondrial and ribosomal genes
remove <- c(grep(Markers_filt$gene, pattern = "^MT-*"), grep(Markers_filt$gene, pattern = "^RPS*"), grep(Markers_filt$gene, pattern = "^RPL*"))
remove
Markers_filt <- Markers_filt[-remove,]

#Select for low % pct2 
Markers_filt <- Markers_filt[(Markers_filt$pct.2 < 0.6),]

#Reorder by high pct.1 and then cell type
Markers_filt <- Markers_filt[order(Markers_filt[,"pct.1"], decreasing = T),] 
Markers_filt <- Markers_filt[order(Markers_filt[,"cluster"]),]
View(Markers_filt)
dim(Markers_filt)
#write.table(Markers_filt, sep = "\t", file = "220425_PS_high_deep_ClusterMarkers_Filtered.txt")
#Markers_filt <- read.table(sep = "\t", file = "PS_deep_glia_ClusterMarkers_Filtered.txt")

DoHeatmap(PS_high_deep, features = c(Markers_filt$gene), slot = "data", draw.lines = T, group.by = "CellType", raster = F) + scale_fill_gradient(low="grey91", high="blue")
plt <- subset(PS_high_deep, idents = c("DenseBulbousAstrocyte", "DenseSmoothAstrocyte"))
DoHeatmap(plt, features = c(Markers_filt$gene), slot = "data", draw.lines = T, group.by = "CellType", raster = F) + scale_fill_gradient(low="grey91", high="blue")


#Identify most specific genes----------------
genes <- Markers_filt$gene
for (i in seq(from = 1, to = length(genes), by = 8)){
  print(VlnPlot(PS_high_deep, group.by = "CellType", same.y.lims = T, features = genes[i:(i+8)]))
}

for (i in seq(from = 1, to = length(genes), by = 8)){
  print(DotPlot(PS_high_deep, group.by = "CellType", features = genes[i:(i+8)])) + theme(axis.text.x = element_text(angle = 45))
}

#Plot on glia from Bhaduri et. al.-------------
#Isolate glial cells------
load(file = "Neocortex_v3.RData")

#Plot cell type markers
DimPlot(Neocortex, label = T)
FeaturePlot(Neocortex, features = "NEUROD2")
FeaturePlot(Neocortex, features = "DLX2")
FeaturePlot(Neocortex, features = "AIF1")
FeaturePlot(Neocortex, features = "FN1")
FeaturePlot(Neocortex, features = "MKI67")
FeaturePlot(Neocortex, features = "PDGFRA")
FeaturePlot(Neocortex, features = "MBP")
FeaturePlot(Neocortex, features = "SOX9")
FeaturePlot(Neocortex, features = "GFAP")
FeaturePlot(Neocortex, features = "SOX2")
FeaturePlot(Neocortex, features = "CRYAB")
FeaturePlot(Neocortex, features = "HOPX")
FeaturePlot(Neocortex, features = "GJA1")
DimPlot(Neocortex, label = T)

#Select glial clusters
Glia <- c(46, 35, 32, 41, 36, 12, 7, 23, 3, 16, 21, 22, 18, 11, 24, 39, 34)
Neocortex_glia1 <- subset(Neocortex, idents = Glia)
DimPlot(Neocortex_glia1, label = T)
Neocortex_glia1

#Remove low quality samples
Idents(Neocortex_glia1) <- "Individual"
Neocortex_glia1 <- subset(Neocortex_glia1, idents = c("GW14", "GW17"), invert = T)
DimPlot(Neocortex_glia1, label = T)
Neocortex_glia1

#Re-process
Neocortex_glia1 <- NormalizeData(object = Neocortex_glia1)
Neocortex_glia1 <- FindVariableFeatures(object = Neocortex_glia1)
Neocortex_glia1 <- ScaleData(object = Neocortex_glia1, vars.to.regress = "Individual")
Neocortex_glia1 <- RunPCA(object = Neocortex_glia1)
ElbowPlot(Neocortex_glia1)
Neocortex_glia1 <- FindNeighbors(object = Neocortex_glia1)
Neocortex_glia1 <- FindClusters(object = Neocortex_glia1)
Neocortex_glia1 <- RunUMAP(object = Neocortex_glia1, dims = 1:20)
DimPlot(object = Neocortex_glia1, reduction = "umap", label = T)

DimPlot(Neocortex_glia1, group.by = "Individual")

FeaturePlot(Neocortex_glia1, features = "NEUROD2")
FeaturePlot(Neocortex_glia1, features = "DLX2")
FeaturePlot(Neocortex_glia1, features = "AIF1")
FeaturePlot(Neocortex_glia1, features = "FN1")
FeaturePlot(Neocortex_glia1, features = "MKI67")
FeaturePlot(Neocortex_glia1, features = "PDGFRA")
FeaturePlot(Neocortex_glia1, features = "MBP")
FeaturePlot(Neocortex_glia1, features = "SOX9")
FeaturePlot(Neocortex_glia1, features = "GFAP")
FeaturePlot(Neocortex_glia1, features = "SOX2")
FeaturePlot(Neocortex_glia1, features = "CRYAB")
FeaturePlot(Neocortex_glia1, features = "HOPX")
FeaturePlot(Neocortex_glia1, features = "GJA1")

DimPlot(object = Neocortex_glia1, reduction = "tsne", label = T)

#Remove non-glial clusters
Remove <- c(11, 15, 0, 3, 17, 9)

Neocortex_glia2 <- subset(Neocortex_glia1, idents = Remove, invert = T)
Neocortex_glia2

#Re-process
Neocortex_glia2 <- NormalizeData(object = Neocortex_glia2)
Neocortex_glia2 <- FindVariableFeatures(object = Neocortex_glia2)
Neocortex_glia2 <- ScaleData(object = Neocortex_glia2, vars.to.regress = "Individual")
Neocortex_glia2 <- RunPCA(object = Neocortex_glia2)
ElbowPlot(Neocortex_glia2)
Neocortex_glia2 <- FindNeighbors(object = Neocortex_glia2)
Neocortex_glia2 <- FindClusters(object = Neocortex_glia2)
Neocortex_glia2 <- RunUMAP(object = Neocortex_glia2, dims = 1:10)
DimPlot(object = Neocortex_glia2, reduction = "umap", label = T)

DimPlot(Neocortex_glia2, group.by = "Individual")

FeaturePlot(Neocortex_glia2, features = "NEUROD2") 
FeaturePlot(Neocortex_glia2, features = "DLX2") 
FeaturePlot(Neocortex_glia2, features = "AIF1")
FeaturePlot(Neocortex_glia2, features = "FN1")
FeaturePlot(Neocortex_glia2, features = "MKI67") 
FeaturePlot(Neocortex_glia2, features = "PDGFRA")
FeaturePlot(Neocortex_glia2, features = "MBP")
FeaturePlot(Neocortex_glia2, features = "SOX9")
FeaturePlot(Neocortex_glia2, features = "GFAP")
FeaturePlot(Neocortex_glia2, features = "SOX2")
FeaturePlot(Neocortex_glia2, features = "CRYAB")
FeaturePlot(Neocortex_glia2, features = "HOPX")
FeaturePlot(Neocortex_glia2, features = "GJA1")

DimPlot(Neocortex_glia4, reduction = "tsne")

#Remove  non-glial clusters
Remove <- c(7, 12, 14)

Neocortex_glia3 <- subset(Neocortex_glia2, idents = Remove, invert = T)
Neocortex_glia3

Neocortex_glia3 <- NormalizeData(object = Neocortex_glia3)
Neocortex_glia3 <- FindVariableFeatures(object = Neocortex_glia3)
Neocortex_glia3 <- ScaleData(object = Neocortex_glia3, vars.to.regress = "Individual")
Neocortex_glia3 <- RunPCA(object = Neocortex_glia3)
ElbowPlot(Neocortex_glia3)
Neocortex_glia3 <- FindNeighbors(object = Neocortex_glia3)
Neocortex_glia3 <- FindClusters(object = Neocortex_glia3)
Neocortex_glia3 <- RunUMAP(object = Neocortex_glia3, dims = 1:10)
DimPlot(object = Neocortex_glia3, reduction = "umap", label = T)

DimPlot(Neocortex_glia3, group.by = "Individual")

FeaturePlot(Neocortex_glia3, features = "NEUROD2") 
FeaturePlot(Neocortex_glia3, features = "DLX2") 
FeaturePlot(Neocortex_glia3, features = "AIF1")
FeaturePlot(Neocortex_glia3, features = "FN1")
FeaturePlot(Neocortex_glia3, features = "MKI67") 
FeaturePlot(Neocortex_glia3, features = "PDGFRA")
FeaturePlot(Neocortex_glia3, features = "MBP")
FeaturePlot(Neocortex_glia3, features = "SOX9")
FeaturePlot(Neocortex_glia3, features = "GFAP")
FeaturePlot(Neocortex_glia3, features = "SOX2")
FeaturePlot(Neocortex_glia3, features = "CRYAB")
FeaturePlot(Neocortex_glia3, features = "HOPX")
FeaturePlot(Neocortex_glia3, features = "GJA1")
FeaturePlot(Neocortex_glia3, features = "FBLN2")
FeaturePlot(Neocortex_glia3, features = "NOG")

#Remove low quality clusters and neurons 
Remove <- c(0, 12, 14) 
Neocortex_glia4 <- subset(Neocortex_glia3, idents = Remove, invert = T)

#Final object including radial glia, astrocytes, and oligodendrocyte lineage
Neocortex_glia4
Neocortex_glia4 <- NormalizeData(object = Neocortex_glia4)
Neocortex_glia4 <- FindVariableFeatures(object = Neocortex_glia4)
Neocortex_glia4 <- ScaleData(object = Neocortex_glia4, vars.to.regress = "Individual")
Neocortex_glia4 <- RunPCA(object = Neocortex_glia4)
ElbowPlot(Neocortex_glia4)
Neocortex_glia4 <- FindNeighbors(object = Neocortex_glia4)
Neocortex_glia4 <- FindClusters(object = Neocortex_glia4, resolution = 0.5)
Neocortex_glia4 <- RunUMAP(object = Neocortex_glia4, dims = 1:10)
DimPlot(object = Neocortex_glia4, reduction = "umap", label = T)

#Assign cell types
Idents(Neocortex_glia4) <- "seurat_clusters"
DimPlot(Neocortex_glia4, label = T)
RG_c <- WhichCells(Neocortex_glia4, idents = c(0, 1, 2, 6, 8))
Astros_c <- WhichCells(Neocortex_glia4, idents = c(7, 5))
OPCs_c <- WhichCells(Neocortex_glia4, idents = c(3, 4))
Oligos_c <- WhichCells(Neocortex_glia4, idents = 9)

RG <- rep(x= "RadialGlia", times = length(RG_c))
Astros <- rep(x= "Astrocytes", times = length(Astros_c))
OPCs <- rep(x= "OPCs", times = length(OPCs_c))
Oligos <- rep(x= "Oligodendrocytes", times = length(Oligos_c))

names(RG) <- RG_c
names(Astros) <- Astros_c
names(OPCs) <- OPCs_c
names(Oligos) <- Oligos_c

CellType <- c(RG, Astros, OPCs, Oligos)
head(CellType, n=100)
Neocortex_glia4$CellType <- CellType

DimPlot(Neocortex_glia4, group.by = "CellType")

#Plot candidate markers on Neocortex_glia4--------
for (i in seq(from = 1, to = length(genes), by = 8)){
  print(DotPlot(Neocortex_glia4, group.by = "CellType", features = genes[i:(i+8)]))
}

#Plot on adult brain from Velmeshev et. al. --------------------------
#Remove cells from ASD patients-----
Adult <- readRDS(file = "autism.rds")
Idents(Adult) <- "diagnosis"
Adult_ctrl <- subset(Adult, idents = "Control")
Adult_ctrl

#Re-process control cells
Adult_ctrl <- NormalizeData(object = Adult_ctrl)
Adult_ctrl <- FindVariableFeatures(object = Adult_ctrl)
Adult_ctrl <- ScaleData(object = Adult_ctrl, vars.to.regress = "individual")
Adult_ctrl <- RunPCA(object = Adult_ctrl, verbose = T)
ElbowPlot(Adult_ctrl)
Adult_ctrl <- FindNeighbors(object = Adult_ctrl)
Adult_ctrl <- FindClusters(object = Adult_ctrl)
Adult_ctrl <- RunUMAP(object = Adult_ctrl, dims = 1:10)
Adult_ctrl <- RunTSNE(object = Adult_ctrl, dims = 1:10)
DimPlot(object = Adult_ctrl, reduction = "umap", label = T)
DimPlot(object = Adult_ctrl, reduction = "tsne", label = T)

#Plot filtered markers on adult data-----------
for (i in seq(from = 1, to = length(genes), by = 8)){
  print(FeaturePlot(Adult_ctrl, reduction = "tsne", features = genes[i:(i+8)]))
}

#Plot on adult brain from Bakken et. al.--------------------------
Adult_bakken <- readRDS("allenM1_10x.rds")
Idents(Adult_bakken) <- "BroadCellType"
Adult_bakken$BroadCellType <- factor(Adult_bakken$BroadCellType, levels = c("Astro L1 FGFR3 SERPINI2", "Astro L1-6 FGFR3 PLCG1", "Astro L1-6 FGFR3 AQP1", "Oligodendrocytes", "Oligo. Progenitors",  "Excitatory Neurons", "Inhibitory Neurons", "Microglia", "Endothelial", "Other"))

#Plot filtered markers on Adult data from Bakken et. al.
for (i in seq(from = 1, to = length(genes), by = 8)){
  print(DotPlot(Adult_bakken, group.by = "BroadCellType", features = genes[i:(i+8)]))
}
