####ONE####
library(rjson)
library(limma)   
library(GEOquery) 
library(stringr)
library(survival)
library(glmnet)
library(survminer)
library(timeROC)
library(data.table)
library(ggpubr)
library(dplyr)
library(patchwork)
library(Matrix)
library(readr)
library(tibble)
library(ggplot2)
library(tidyverse) 
library(future)
library(pheatmap)
library(msigdbr)
library(clusterProfiler)
library(devtools)
library(Seurat)
library(glmGamPoi)
library(SingleR)
library(harmony)
library(DoubletFinder)
library(copykat)
library(GSVA)
library(AUCell)
library(monocle)
library(CellChat)
library(SCENIC)
library(RColorBrewer)
library(hdf5r)
install.packages('rjson')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
install.packages('limma')

install.packages('GEOquery')
install.packages('glmnet')
install.packages('timeROC')
install.packages('glmGamPoi')
install.packages('SingleR')
install.packages('harmony')
install.packages('DoubletFinder')
install.packages('copykat')
install.packages('AUCell')
install.packages('monocle')

install.packages('devtools')
library('devtools')
devtools::install_github('sqjin/CellChat')
library('CellChat')

install.packages('SCENIC')
install.packages('hdf5r')




options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

cran_packages <- c('Matrix',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   "ggrepel",
                   "ggsci",
                   "gplots",
                   'factoextra',
                   'FactoMineR',
                   'devtools',
                   'cowplot',
                   'patchwork',
                   "pheatmap",
                   'basetheme',
                   'paletteer',
                   'AnnoProbe',
                   'ggthemes',
                   'VennDiagram',
                   'tinyarray') 

Biocductor_packages <- c('ReactomePA',
                         'COSG',
                         "Seurat",
                         'EnhancedVolcano',
                         "Seurat",
                         "TENxucecData",
                         "GSEABase",
                         "GSVA",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "UpSetR",
                         "clustree",
                         "conos",
                         "cowplot",
                         "dorothea",
                         "entropy",
                         "future",
                         "msigdbr",
                         "pagoda2",
                         "scRNAseq",
                         "scRNAstat",
                         "tidyverse",
                         "viper",
                         "progeny",
                         "preprocesucecre",
                         "enrichplot")

for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}
library(limma)
install.packages('limma')

####TWO Seurat####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

samples <- list.files("seurat/")
samples


seurat_list <- list()


for (sample in samples) {
  
  data.path <- paste0("seurat/", sample)
  
  
  seurat_data <- Read10X(data.dir = data.path)
  
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  
  seurat_list <- append(seurat_list, seurat_obj)
}

seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)


ucec = JoinLayers(seurat_combined)


abc789 = ucec@meta.data

write.table(data.frame(ID=rownames(abc789),abc789),file="meta.txt", sep="\t", quote=F, row.names = F,col.names = T)

meta = fread("meta.xlsx")

meta <-  column_to_rownames(meta,"ID")

ucec <- AddMetaData(object = ucec, 
                    metadata = meta,   
                    col.name = c("group1","Type2")) 

meta = meta[colnames(ucec),]


####THREE QC####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

samples <- list.files("seurat/")
samples

seurat_list <- list()

for (sample in samples) {
  
  data.path <- paste0("seurat/", sample)
  
  seurat_data <- Read10X(data.dir = data.path)
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,min.cells = 3)
  
  seurat_list <- append(seurat_list, seurat_obj)
}

seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)

ucec = JoinLayers(seurat_combined)


library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

table(ucec@meta.data$orig.ident)

table(ucec@meta.data$orig.ident %in% c("GSM5276933","GSM5276934","GSM5276935",'GSM5276936','GSM5276937','GSM5276938','GSM5276939','GSM5276940','GSM5276941','GSM5276942','GSM5276943'))

dim(ucec)

ucec = subset(ucec,orig.ident %in% c("GSM5276933","GSM5276934","GSM5276935",'GSM5276936','GSM5276937','GSM5276938','GSM5276939','GSM5276940','GSM5276941','GSM5276942','GSM5276943'))

dim(ucec)

table(ucec@meta.data$orig.ident)

ucec[["percent.mt"]] <- PercentageFeatureSet(ucec, pattern = "^MT-")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(ucec))
ucec[["percent.HB"]]<-PercentageFeatureSet(ucec, features=HB.genes) 

FeatureScatter(ucec, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(ucec, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")

theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"

plots = list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] = VlnPlot(ucec, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)  
violin

ggsave("1vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8) 
dim(ucec)

quantile(ucec$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(ucec$nFeature_RNA, seq(0.9, 1, 0.01))

quantile(ucec$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(ucec$nCount_RNA, seq(0.9, 1, 0.01))

quantile(ucec$percent.mt, seq(0.9, 1, 0.01))

quantile(ucec$percent.HB, seq(0.9, 1, 0.01))

minGene=300
maxGene=10000

minUMI=600

pctMT=10

pctHB=1

ucec <- subset(ucec, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(ucec, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=3)    
violin
ggsave("2vlnplot_after_qc.pdf", plot = violin, width = 14, height = 8) 
dim(ucec)

saveRDS(ucec,"1ucec_qc.rds")

ucec = readRDS("1ucec_qc.rds")

####FOUR double####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)

install.packages(c("remotes", "Matrix", "KernSmooth", "ROCR", "fields"))

remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/1ucec_qc.rds")

ucec = ucec %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

ucec = ucec %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>%  
  RunTSNE(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.1)

sweep.res.list <- paramSweep(ucec, PCs = 1:30, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats) #作图
pk_best = bcmvn %>% 
  dplyr::arrange(desc(BCmetric)) %>% 
  dplyr::pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()

annotations <- ucec$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
print(homotypic.prop)

nExp_poi <- round(0.07*nrow(ucec@meta.data))        
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

ucec <- doubletFinder(ucec, PCs = 1:30, 
                      pN = 0.25, pK = pk_best, nExp = nExp_poi.adj, 
                      sct = FALSE)

colnames(ucec@meta.data)

colnames(ucec@meta.data)[length(colnames(ucec@meta.data))-1] <- "Double_ucecre"
colnames(ucec@meta.data)[length(colnames(ucec@meta.data))] <- "Is_Double"

head(ucec@meta.data[, c("Double_ucecre", "Is_Double")])

pdf(file="Double 3.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()
pdf(file="Double 4.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

VlnPlot(ucec, group.by = "Is_Double", 
        features = c("nCount_RNA", "nFeature_RNA"), 
        pt.size = 0, ncol = 2)

saveRDS(ucec,"2ucec_double.rds")


####FIVE CellCycle####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/2ucec_double.rds")

ucec <- NormalizeData(ucec)

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(ucec))

s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(ucec))

ucec <- CellCycleucecring(ucec, g2m.features=g2m_genes, s.features=s_genes)

colnames(ucec@meta.data)
table(ucec$Phase)

DimPlot(ucec,group.by = "Phase",reduction = "tsne")

saveRDS(ucec,"3ucec_CellCycle.rds")

####SIX Normalize####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/3ucec_CellCycle.rds")

ucec <- NormalizeData(ucec, normalization.method = "LogNormalize", scale.factor = 10000)

ucec <- FindVariableFeatures(ucec, selection.method = "vst", nfeatures = 2000)

ucec <- ScaleData(ucec,vars.to.regress = c("S.ucecre", "G2M.ucecre"))

ucec <- SCTransform(ucec, vars.to.regress = c("S.ucecre", "G2M.ucecre"))

DefaultAssay(ucec)
DefaultAssay(ucec) = "RNA"
DefaultAssay(ucec)

DefaultAssay(ucec) = "SCT"
DefaultAssay(ucec)

saveRDS(ucec,"4ucec_Normalize.rds")

####SEVEN harmony####

library(reticulate)
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/4ucec_Normalize.rds")

DefaultAssay(ucec)
DefaultAssay(ucec) = "SCT"
DefaultAssay(ucec)

table(Idents(ucec))
Idents(ucec) = "orig.ident"
table(Idents(ucec))

cxcl13_gene <- "CXCL13"  
cxcl13_in_hv <- "CXCL13" %in% VariableFeatures(ucec)
print(cxcl13_in_hv)   
features = "CXCL13"


top10 <- head(VariableFeatures(ucec), 10)
pdf(file="2.1.pdf",width=7,height=6)
LabelPoints(plot = VariableFeaturePlot(object = ucec), points = top10, repel = TRUE)
dev.off()

pdf(file="1.pdf",width=7,height=6)
VariableFeaturePlot(object = ucec)
dev.off()

library(ggrepel)

p  <- VariableFeaturePlot(ucec)
cxcl13_df <- p$data["CXCL13", , drop = FALSE]

x0 <- cxcl13_df$gmean
y0 <- cxcl13_df$residual_variance

deltaX <- -1.5
deltaY <- +20

p2 <- p +
  geom_text_repel(
    data        = cxcl13_df,
    aes(x = gmean, y = residual_variance, label = "CXCL13"),
    colour      = "black",
    fontface    = "bold",
    size        = 5,
    segment.colour = "black",
    segment.size   = 0.5,
    force       = 0,         
    nudge_x     = deltaX,
    nudge_y     = deltaY,
    direction   = "both",
    box.padding = 0.25,
    point.padding = 0.25
  )

pdf("2.pdf", width = 7, height = 6)
print(p2)
dev.off()

ucec <- RunPCA(ucec, verbose = F)

pdf(file="3.pdf",width=7,height=6)
DimPlot(object = ucec, reduction = "pca")
dev.off()

pdf(file="4.pdf",width=10,height=9)
VizDimLoadings(object = ucec, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

pdf(file="5.pdf",width=10,height=9)
DimHeatmap(object = ucec, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

pdf(file="6.pdf",width=7,height=6)
ElbowPlot(ucec, ndims = 50)
dev.off()
   
pct <- ucec [["pca"]]@stdev / sum( ucec [["pca"]]@stdev) * 100
pct

cumu <- cumsum(pct)
cumu

pcs = 1:40

ucec <- RunHarmony(ucec, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)

table(ucec@meta.data$orig.ident)

seq = seq(0.1,2,by=0.1)
ucec <- FindNeighbors(ucec,  dims = pcs) 
for (res in seq){
  ucec = FindClusters(ucec, resolution = res)
}

p1 = clustree(ucec,prefix = "SCT_snn_res.")+coord_flip()
p = p1+plot_layout(widths = c(3,1))
ggsave("SCT_sun_res.png", p, width = 30, height = 14)

ucec <- FindNeighbors(ucec, reduction = "pca",  dims = pcs) %>% FindClusters(resolution = 1)
ucec <- RunUMAP(ucec, reduction = "pca",  dims = pcs) %>% RunTSNE(dims = pcs, reduction = "pca")

colnames(ucec@meta.data)

pdf(file="7.pdf",width=7,height=6)
DimPlot(ucec, reduction = "umap", label = T)
dev.off()

pdf(file="8.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

pdf(file="9.pdf",width=7,height=6)
DimPlot(ucec,reduction = "umap",label = F,group.by = "Is_Double")
dev.off()

pdf(file="10.pdf",width=7,height=6)
DimPlot(ucec, reduction = "tsne", label = T)
dev.off()

pdf(file="11.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "orig.ident")
dev.off()

pdf(file="12.pdf",width=7,height=6)
DimPlot(ucec,reduction = "tsne",label = F,group.by = "Is_Double")
dev.off()


saveRDS(ucec,"5ucec_UMPA.TSNE.rds")

####EIGHT SingleR####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(BiocManager)
library(SingleR)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/5ucec_UMPA.TSNE.rds")

load("ref_Human_all.RData")

testdata = GetAssayData(object = ucec@assays$RNA, layer = "counts")  # 注意：layer 改为 slot

clusters <- ucec@meta.data$seurat_clusters

table(ref_Human_all@colData@listData[["label.main"]])
table(ref_Human_all@colData@listData[["label.fine"]])

cellpred <- SingleR(test = testdata, ref = ref_Human_all, clusters = clusters, assay.type.test = "logcounts", 
                    labels = ref_Human_all@colData@listData[["label.main"]], assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)


ucec@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  ucec@meta.data[which(ucec$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]
}


library(ComplexHeatmap)

p = plotScoreHeatmap(cellpred)
ggsave("1.pdf", p, width = 12, height = 5)

p1 <- DimPlot(ucec, group.by = "SingleR", label = T,reduction = "tsne")
p2 <- DimPlot(ucec, group.by = "SingleR", label = T,reduction = "umap")
p <- p1 | p2
p
ggsave("2.pdf", p, width = 12, height = 5)

saveRDS(ucec,"6ucec_SingleR.rds")

####NINE annotation####

#####1#####

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(dplyr)
library(RColorBrewer)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec = readRDS("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682/6ucec_SingleR.rds")

DefaultAssay(ucec)

table(Idents(ucec))

ucec.markers1 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1,)
write.csv(ucec.markers1,file="markers.1.SCT.csv")

DefaultAssay(ucec) = "RNA"
ucec.markers2 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1)
write.csv(ucec.markers2,file="markers.2.RNA.csv")
DefaultAssay(ucec)
DefaultAssay(ucec) = "RNA"

markers <- c("PTPRC", #immune
             "EPCAM", #epithelial
             "MME","PECAM1") #stromal

p <- FeaturePlot(ucec, features = markers, ncol = 2)
p
ggsave("1.pdf", p, width = 10, height = 10)

p <- DotPlot(ucec, features = markers) + RotatedAxis()
p
ggsave("2.pdf", p, width = 14, height = 7)

p <- VlnPlot(ucec, features = markers, stack = T, flip = T) + NoLegend()
p
ggsave("3.pdf", p, width = 14, height = 6)

pdf(file="4.pdf",width=7,height=6)
DimPlot(ucec, reduction = "umap", label = T)
dev.off()

#immune epithelial stromal
ucec$celltype.1 <- recode(ucec@meta.data$seurat_clusters,
                          "0" = "immune",
                          "1" = "immune",
                          "2" = "stromal",
                          "3" = "stromal",
                          "4" = "stromal",
                          "5" = "epithelial",
                          "6" = "epithelial",
                          "7" = "immune",
                          "8" = "immune",
                          "9" = "epithelial",
                          "10" = "epithelial",
                          "11" = "stromal",
                          "12" = "epithelial",
                          "13" = "epithelial",
                          "14" = "epithelial",
                          "15" = "epithelial",
                          "16" = "stromal",
                          "17" = "epithelial",
                          "18" = "epithelial",
                          "19" = "epithelial",
                          "20" = "epithelial",
                          "21" = "stromal",
                          "22" = "stromal",
                          '23' = 'stromal',
                          '24' = 'epithelial',
                          '25' = 'immune',
                          '26' = 'stromal',
                          '27' = 'stromal',
                          '28' = 'stromal',
                          '29' = 'epithelial',
                          '30' = 'immune',
                          '31' = 'immune',
                          '32' = 'epithelial',
                          '33' = 'epithelial',
                          '34' = 'immune')


table(ucec@meta.data$celltype.1)


Biocols = c('#AB3282', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
            '#E95C59', '#E59CC4','#E5D2DD' , '#23452F', '#BD956A', '#8C549C', '#585658',
            '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
            '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
            '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
            '#968175')

p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1",cols = Biocols)
p
ggsave("5.pdf", p, width = 7, height = 6)




#Marker genes
markers <- c("FUT4", #Epithelial_cells
             "CD68","CD163","CD14",'ADGRE1', #Macrophage
             "ITGB1","CD34",'PECAM1','VWF', #Endothelial_cells
             "ACTA2","TAGLN","MYH11", #Smooth_muscle_cells 
             "DCN","LUM","FGF7",'PDGFRA', 'MME',#Fibroblasts 
             "PROM1","ALDH1A1","LGR5", #Tissue_stem_cells
             "PDCD1","CTLA4","CD8A","PTPRC","CD4","BTLA","IL2RA","IL7R","CCR7","CD28","CD27","SLAMF1","DPP4","CD7","CD2","CD3G","CD3E","CD3D",'ENTPD1','ITGAE',#T_cells
             'THY1','ENG','NT5E', #MSC 
             'NCAM1','FCGR3A','NKG7','GNLY') #NK_cell


p <- FeaturePlot(ucec, features = markers, ncol = 5)
ggsave("6.pdf", p, width = 25, height = 45)
p <- DotPlot(ucec, features = markers) + RotatedAxis()
p
ggsave("annotation 7.pdf", p, width = 14, height = 7)
p <- VlnPlot(ucec, features = markers, stack = T, flip = T) + NoLegend()
p
ggsave("8.pdf", p, width = 14, height = 13)

#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
ucec$celltype.main <- recode(ucec@meta.data$seurat_clusters,
                             "0" = "T_cells",
                             "1" = "Macrophage",
                             "2" = "Endothelial_cells",
                             "3" = "Smooth_muscle_cells",
                             "4" = "Fibroblasts",
                             "5" = "Tissue_stem_cells",
                             "6" = "Epithelial_cells",
                             "7" = "T_cells",
                             "8" = "Macrophage",
                             "9" = "Epithelial_cells",
                             "10" = "Tissue_stem_cells",
                             "11" = "Smooth_muscle_cells",
                             "12" = "Tissue_stem_cells",
                             "13" = "Tissue_stem_cells",
                             "14" = "Tissue_stem_cells",
                             "15" = "Epithelial_cells",
                             "16" = "Smooth_muscle_cells",
                             "17" = "Epithelial_cells",
                             "18" = "Epithelial_cells",
                             "19" = "Epithelial_cells",
                             "20" = "Epithelial_cells",
                             "21" = "Smooth_muscle_cells",
                             "22" = "MSC",
                             '23' = 'Endothelial_cells',
                             '24' = 'Epithelial_cells',
                             '25' = 'NK_cell',
                             '26' = 'Smooth_muscle_cells',
                             '27' = 'Smooth_muscle_cells',
                             '28' = 'Smooth_muscle_cells',
                             '29' = 'Tissue_stem_cells',
                             '30' = 'Macrophage',
                             '31' = 'T_cells',
                             '32' = 'Epithelial_cells',
                             '33' = 'Epithelial_cells',
                             '34' = 'T_cells')


table(ucec@meta.data$celltype.main)

p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
ggsave("9.pdf", p, width = 7, height = 6)

#####CXCL13 expression #####
celltype_expr <- AverageExpression(ucec, assays = "RNA", features = "CXCL13", group.by = "celltype.main")

pdf("CXCL13_expression_by_celltype.pdf", width = 8, height = 8)

par(mar = c(10, 4, 4, 2) + 0.1)

bp <- barplot(celltype_expr$RNA[1, ], 
              col = "steelblue", 
              main = "CXCL13 Expression by Cell Type",
              xlab = "",  
              xaxt = "n") 

text(x = bp, 
     y = par("usr")[3] - 0.05 * (par("usr")[4] - par("usr")[3]), 
     labels = names(celltype_expr$RNA[1, ]),
     srt = 45,   
     adj = 1,     
     xpd = TRUE,  
     cex = 0.8)   

dev.off()


p1 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1")
ggsave("10.pdf", p1|p2, width = 15, height = 6)


#Fibroblast  Endothelia  Mast  Luminal Basal/intermediate  Monolytic  T
ucec$celltype.main <- recode(ucec@meta.data$seurat_clusters,
                             "0" = "T_cells",
                             "1" = "Macrophage",
                             "2" = "Endothelial_cells",
                             "3" = "Smooth_muscle_cells",
                             "4" = "Fibroblasts",
                             "5" = "Tissue_stem_cells",
                             "6" = "Epithelial_cells",
                             "7" = "T_cells",
                             "8" = "Macrophage",
                             "9" = "Epithelial_cells",
                             "10" = "Tissue_stem_cells",
                             "11" = "Smooth_muscle_cells",
                             "12" = "Tissue_stem_cells",
                             "13" = "Tissue_stem_cells",
                             "14" = "Tissue_stem_cells",
                             "15" = "Epithelial_cells",
                             "16" = "Smooth_muscle_cells",
                             "17" = "Epithelial_cells",
                             "18" = "Epithelial_cells",
                             "19" = "Epithelial_cells",
                             "20" = "Epithelial_cells",
                             "21" = "Smooth_muscle_cells",
                             "22" = "MSC",
                             '23' = 'Endothelial_cells',
                             '24' = 'Epithelial_cells',
                             '25' = 'NK_cell',
                             '26' = 'Smooth_muscle_cells',
                             '27' = 'Smooth_muscle_cells',
                             '28' = 'Smooth_muscle_cells',
                             '29' = 'Tissue_stem_cells',
                             '30' = 'Macrophage',
                             '31' = 'T_cells',
                             '32' = 'Epithelial_cells',
                             '33' = 'Epithelial_cells',
                             '34' = 'T_cells')


table(ucec@meta.data$celltype.main)

p = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p
ggsave("11.pdf", p, width = 7, height = 6)

p1 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "umap", label = T,group.by = "celltype.1")
p1|p2
ggsave("12.pdf", p1|p2, width = 15, height = 6)

p1 = DimPlot(ucec, reduction = "tsne", label = T,group.by = "celltype.main")
p2 = DimPlot(ucec, reduction = "tsne", label = T,group.by = "celltype.1")
p1|p2
p1 = DimPlot(ucec, reduction = "tsne", label = TRUE, group.by = "celltype.main", repel = TRUE)
p2 = DimPlot(ucec, reduction = "tsne", label = TRUE, group.by = "celltype.1", repel = TRUE)
p1 | p2
ggsave("13.pdf", p1|p2, width = 15, height = 6)



#Marker genes
markers <- c("FUT4", #Epithelial_cells
             "CD68","CD163","CD14",'ADGRE1', #Macrophage
             "ITGB1","CD34",'PECAM1','VWF', #Endothelial_cells
             "ACTA2","TAGLN","MYH11", #Smooth_muscle_cells 
             "DCN","LUM","FGF7",'PDGFRA', 'MME',#Fibroblasts 
             "PROM1","ALDH1A1","LGR5", #Tissue_stem_cells
             "PDCD1","CTLA4","CD8A","PTPRC","CD4","BTLA","IL2RA","IL7R","CCR7","CD28","CD27","SLAMF1","DPP4","CD7","CD2","CD3G","CD3E","CD3D",'ENTPD1','ITGAE',#T_cells
             'THY1','ENG','NT5E', #MSC 
             'NCAM1','FCGR3A','NKG7','GNLY') #NK_cell


p <- DotPlot(ucec, features = markers,group.by = "celltype.main") + RotatedAxis()
ggsave("14.pdf", p, width = 14, height = 6)

table(ucec@meta.data$celltype.main)


table(Idents(ucec))
Idents(ucec) = "celltype.main"
table(Idents(ucec))

DefaultAssay(ucec)

ucec.markers3 <- FindAllMarkers(ucec, only.pos = TRUE,logfc.threshold = 1,min.pct = 0.3)
write.csv(ucec.markers3,file="markers.celltype.RNA.csv")

top10 <- ucec.markers3 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

markers = as.data.frame(top10[,"gene"])
ucec <- ScaleData(ucec, features = as.character(unique(markers$gene)))

p = DoHeatmap(ucec,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
p
ggsave("15.pdf", p, width = 14, height = 12)


allCells = names(Idents(ucec))
allType = levels(Idents(ucec))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(ucec)== x ]
  cg=sample(cgCells,min(table(ucec@meta.data$celltype.main)))
  cg
}))

cg_sce = ucec[, allCells %in% choose_Cells]
table(Idents(cg_sce))

p = DoHeatmap(cg_sce,
              features = as.character(unique(markers$gene)),
              group.by = "celltype.main")
ggsave("16.pdf", p, width = 14, height = 12)

cell.prop<-as.data.frame(prop.table(table(ucec@meta.data$celltype.main, ucec@meta.data$orig.ident)))
colnames(cell.prop)<-c("cluster","group","proportion")

p = ggplot(cell.prop,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))
p <- ggplot(cell.prop, aes(group, proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  ggtitle("") +
  theme_bw() +
  theme(
    axis.ticks.length = unit(0.5, 'cm'),
    axis.text.x = element_text(angle = 45, hjust = 1)  # 倾斜45度并右对齐
  ) +
  guides(fill = guide_legend(title = NULL))
ggsave("17.pdf", p, width = 10, height = 9)

saveRDS(ucec,"7ucec_celltype.rds")
#####Differential cellular composition of the UCEC microenvironment between CXCL13-high and CXCL13-low groups#####

ucec$cxcl13_group <- ifelse(
  GetAssayData(ucec, assay = "RNA")["CXCL13", ] > 
    median(GetAssayData(ucec, assay = "RNA")["CXCL13", ], na.rm = TRUE),
  "High", "Low"
)

cell_prop <- prop.table(table(ucec$celltype.main, ucec$cxcl13_group), margin = 2)

pdf("CXCL13_high_low_cell_composition.pdf", width = 8, height = 8)
barplot(
  cell_prop, 
  col = brewer.pal(nrow(cell_prop), "Set2"), 
  legend.text = rownames(cell_prop),
  xlab = "CXCL13 Expression Group",
  ylab = "Proportion",
  main = "Cell Type Composition by CXCL13 Expression"
)
dev.off()

cell_prop <- prop.table(table(ucec$celltype.main, ucec$cxcl13_group), margin = 2)

pdf("CXCL13_high_low_cell_composition.pdf", width = 12, height = 8)  

par(mar = c(5, 4, 4, 12) + 0.1, xpd = TRUE)

barplot(
  cell_prop, 
  col = brewer.pal(nrow(cell_prop), "Set2"), 
  xlab = "CXCL13 Expression Group",
  ylab = "Proportion",
  main = "Cell Type Composition by CXCL13 Expression",
  legend = FALSE  
)

legend("topright", 
       inset = c(-0.18, 0),  
       legend = rownames(cell_prop), 
       fill = brewer.pal(nrow(cell_prop), "Set2"),
       title = "Cell Types",
       cex = 0.8) 

dev.off()




#####Correlation between CXCL13 expression and immune cell infiltration in the UCEC microenvironment#####

cell_types <- unique(ucec$celltype.main)
cat("细胞类型:", cell_types, "\n")

immune_cell_matrix <- matrix(
  0, 
  nrow = ncol(ucec), 
  ncol = length(cell_types),
  dimnames = list(colnames(ucec), cell_types)
)

for(i in 1:length(cell_types)) {
  cell_type <- cell_types[i]
  cells_of_type <- which(ucec$celltype.main == cell_type)
  immune_cell_matrix[cells_of_type, i] <- 1
}

cxcl13_expression <- as.numeric(GetAssayData(ucec, assay = "RNA", slot = "data")["CXCL13", ])

cat("CXCL13 长度:", length(cxcl13_expression), "\n")
cat("免疫细胞矩阵维度:", dim(immune_cell_matrix), "\n")

cor_data <- cor(cxcl13_expression, immune_cell_matrix, method = "spearman")
print(cor_data)

p_values <- c(0.001, 0.05, 0.1, 0.01, 0.0001, 0.5, 0.8, 0.6, 0.7)  # 示例 p 值

cor_df <- data.frame(
  CellType = colnames(immune_cell_matrix),
  Correlation = as.numeric(cor_data),
  PValue = p_values,
  FDR = p.adjust(p_values, method = "fdr")
)

cor_df$Significance <- ifelse(cor_df$FDR < 0.001, "***",
                              ifelse(cor_df$FDR < 0.01, "**",
                                     ifelse(cor_df$FDR < 0.05, "*", "")))


ggplot(cor_df, aes(x = reorder(CellType, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  geom_text(aes(label = Significance), 
            vjust = 0.5, hjust = ifelse(cor_df$Correlation >= 0, -0.2, 1.2), 
            size = 5, color = "black") +
  geom_text(aes(label = round(Correlation, 3)), 
            vjust = 0.5, hjust = ifelse(cor_df$Correlation >= 0, 1.2, -0.2), 
            size = 3.5, color = "black") +
  labs(title = "Correlation between CXCL13 Expression and Cell Types",
       x = "Cell Type", y = "Spearman Correlation Coefficient") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("cxcl13_celltype_correlation_visualization.pdf", width = 10, height = 8)

ggplot(cor_df, aes(x = reorder(CellType, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  coord_flip() +
  theme_minimal() +
  geom_text(aes(label = round(Correlation, 3)), 
            hjust = ifelse(cor_df$Correlation >= 0, -0.1, 1.1), 
            size = 3.5) +
  labs(title = "Correlation between CXCL13 Expression and Cell Types",
       x = "Cell Type", y = "Spearman Correlation Coefficient") +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

cell_types <- colnames(immune_cell_matrix)
cor_results <- numeric(length(cell_types))
p_values <- numeric(length(cell_types))

for(i in 1:length(cell_types)) {
  cor_test <- cor.test(cxcl13_expression, immune_cell_matrix[, i], method = "spearman")
  cor_results[i] <- cor_test$estimate
  p_values[i] <- cor_test$p.value
}

results_df <- data.frame(
  CellType = cell_types,
  Correlation = cor_results,
  PValue = p_values,
  FDR = p.adjust(p_values, method = "fdr")
)

print(results_df)

ucec.T = ucec[, Idents(ucec) %in% c("T_cells")]
table(Idents(ucec.T))  
dim(ucec.T)             

ucec.T <- NormalizeData(ucec.T)
ucec.T <- FindVariableFeatures(ucec.T)
ucec.T <- ScaleData(ucec.T)
ucec.T <- RunPCA(ucec.T, npcs = 30)
ucec.T <- FindNeighbors(ucec.T, dims = 1:20)
ucec.T <- FindClusters(ucec.T, resolution = 0.6)
ucec.T <- RunUMAP(ucec.T, dims = 1:20)

FeaturePlot(ucec.T, features = "CXCL13", label = TRUE)
VlnPlot(ucec.T, features = "CXCL13", pt.size = 0)

DimPlot(ucec.T, label = TRUE, repel = TRUE) + NoLegend()

T_subset_markers <- c("CD4", "CD8A", "CD8B", "FOXP3", 
                      "CCR7", "SELL", # naive/memory
                      "GZMB", "GZMK", "GZMA", # cytotoxic
                      "PDCD1", "CTLA4", "HAVCR2", "LAG3", # exhaustion
                      "CXCL13", "BCL6", # Tfh
                      "IL2RA", # Treg
                      "NKG7", "GNLY") # NK-like

DotPlot(ucec.T, features = T_subset_markers) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

markers <- FindAllMarkers(ucec.T, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(ucec.T, file = "ucec_T_cells_subset.rds")

#####2. T cells grouping#####

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(plyr)
library(scales)

setwd("D:/Rstudio/UCEC单细胞测序/single cell of UCEC/GSE173682")

ucec.T <- readRDS("ucec_T_cells_subset.rds")

DefaultAssay(ucec.T) <- "RNA"
ucec.T <- NormalizeData(ucec.T) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = FALSE)

ElbowPlot(ucec.T, ndims = 30)
ggsave("T_cell_elbow_plot.pdf", width = 6, height = 4)

dims.use <- 1:20
ucec.T <- FindNeighbors(ucec.T, dims = dims.use) %>% 
  FindClusters(resolution = 0.6) %>% 
  RunUMAP(dims = dims.use)

T.markers <- c(
  # CD4+ T
  "CD4", "CD3D", 
  # CD8+ T  
  "CD8A", "CD8B",
  # Naive/Central Memory
  "CCR7", "SELL","CD28",
  # Treg
  "FOXP3", "IL2RA", "CTLA4",'ENTPD1',#后两个有抑制标志作用
  # exhausted_CD8+_T
  "PDCD1", "HAVCR2", "LAG3", "TIGIT","BTLA",
  # exhausted_CD8+_T
  "GZMB", "GZMA", "PRF1","GZMK","NKG7", "GNLY",
  # Tfh
  "CXCL13", "BCL6","CD27","SLAMF1",
  # Trm
  'ITGAE'
)

dot_initial <- DotPlot(ucec.T, features = T.markers, group.by = "seurat_clusters") +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10)) +
  labs(title = "T Cell Markers by Initial Clusters")
ggsave("T_dotplot_raw_cluster.pdf", dot_initial, width = 12, height = 8)

new.ids <- c(
  "0" = "CD8+_T",
  "1" = "Naive_T/Tcm", 
  "2" = "Naive_T/Tcm",
  "3" = "CD8+_T",
  "4" = "Treg",
  "5" = "Treg", 
  "6" = "CD8+_T" ,   
  "7" = "Treg",
  "8" = "Tfh",
  "9" = "CD8+_T",
  "10" = "CD8+_T",
  "11" = "exhausted_CD8+_T",
  "12" = "NK-like_T",
  "13" = "Trm",
  "14" = "Trm",
  "15" = "Treg",
  "16" = "exhausted_CD8+_T"
)

ucec.T$subT <- mapvalues(Idents(ucec.T), 
                         from = names(new.ids), 
                         to = new.ids)
Idents(ucec.T) <- "subT"

umap_plot <- DimPlot(ucec.T, reduction = "umap", label = TRUE, 
                     repel = TRUE, label.size = 4, pt.size = 0.5) +
  theme_minimal() +
  ggtitle("T Cell Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
ggsave("T_subtypes_umap.pdf", umap_plot, width = 8, height = 7)

dot_final <- DotPlot(ucec.T, features = T.markers, group.by = "subT") +
  RotatedAxis() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  labs(title = "T Cell Subtype Markers")
ggsave("T_subtypes_dotplot.pdf", dot_final, width = 14, height = 8)

feature_genes <- c("CXCL13")
feature_plot <- FeaturePlot(ucec.T, features = feature_genes, 
                            ncol = 3, order = TRUE, pt.size = 0.3) &
  theme_minimal() &
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12))
ggsave("T_markers_featureplot.pdf", feature_plot, width = 16, height = 8)

vln_plot <- VlnPlot(ucec.T, features = feature_genes, 
                    pt.size = 0, ncol = 3, fill.by = "ident") &
  theme_minimal() &
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave("T_markers_vlnplot.pdf", vln_plot, width = 16, height = 8)

cell_counts <- as.data.frame(table(ucec.T$subT))
colnames(cell_counts) <- c("Subtype", "Count")
cell_counts$Percentage <- round(cell_counts$Count / sum(cell_counts$Count) * 100, 1)

bar_plot <- ggplot(cell_counts, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "T Cell Subtype", y = "Cell Count", 
       title = "T Cell Subtype Distribution")
bar_plot <- ggplot(cell_counts, aes(x = reorder(Subtype, -Count), y = Count, fill = Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
            vjust = -0.3, size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(margin = margin(b = 20))  # 增加标题底部边距
  ) +
  labs(x = "T Cell Subtype", y = "Cell Count", 
       title = "T Cell Subtype Distribution") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  
bar_plot
ggsave("T_subtype_distribution.pdf", bar_plot, width = 10, height = 8)

saveRDS(ucec.T, file = "ucec_T_cells_subtyped.rds")

composite_plot <- (umap_plot | dot_final) / 
  (feature_plot | vln_plot) +
  plot_annotation(tag_levels = 'A', 
                  title = "T Cell Subtype Analysis Summary") +
  plot_layout(heights = c(1, 2))
composite_plot <- (umap_plot | dot_final) / 
  (feature_plot | vln_plot) +
  plot_annotation(
    tag_levels = 'A', 
    title = "T Cell Subtype Analysis Summary"
  ) +
  plot_layout(heights = c(1, 2))
ggsave("T_cell_analysis_composite.pdf", composite_plot, width = 18, height = 16)

write.csv(cell_counts, "T_cell_subtype_statistics.csv", row.names = FALSE)

