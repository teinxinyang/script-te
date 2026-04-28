library(Seurat)
library(tidyverse)

nsclc.sparse.m <- Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5')
names(nsclc.sparse.m)
cts <-  nsclc.sparse.m$`Gene Expression`

# Read10X_h5 não devolve apenas uma tabela de números. 
# Ela devolve uma Lista Organizada (pense em uma pasta com várias divisórias).
# Vá até a variável nsclc.sparse.m, abra a divisória chamada Gene Expression e 
# salve apenas o que está lá dentro em uma nova variável chamada cts

nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
# pelo menos células que expressem 200 genes diferentes
str(nsclc.seurat.obj)
# ver como ele é organizado por dentro, quais são os compartimentos e que tipo de dados estão lá."
nsclc.seurat.obj
# 1 assay: Significa que você tem 1 modalidade de dado (neste caso, apenas RNA).
# cada uma das 1.575 células que aparecem no seu resumo tem, obrigatoriamente, 200 ou mais genes diferentes detectados nela.

# ================= controle de qualidade ============================================
View(nsclc.seurat.obj@meta.data)
# O nFeature_RNA representa o número de genes diferentes que foram detectados naquela célula específica.
# O nCount_RNA representa o número total de moléculas de RNA (também chamadas de UMIs) detectadas naquela célula, somando todos os genes.
# Na analogia: É o número total de livros físicos na estante. Se você tem 5 exemplares do livro de Farmacologia e 10 de Imunologia, seu nCount é 15, mesmo que você só tenha 2 nFeatures (títulos).



# Se o sequenciador lê uma gota e percebe que uma porcentagem gigante de todo o RNA lá dentro vem de mitocôndrias,
# isso é um atestado de óbito.
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)











