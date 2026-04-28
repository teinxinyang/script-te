import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # sem display gráfico
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white')

print("=" * 50)
print("  Análise scRNA-seq — NSCLC (10x Genomics 5')")
print("=" * 50)

# --- 1. Carregar dados ---
print("\n[1/9] Carregando dados...")
adata = sc.read_10x_mtx('/home/lait/tei/teinxin_tc/vdj_v1_hs_nsclc_multi_5gex_t_b_count_filtered_feature_bc_matrix (1) (1)', var_names='gene_symbols', cache=True)


adata.var_names_make_unique()
print(f"  {adata.n_obs} células, {adata.n_vars} genes")

# --- 2. QC ---
print("\n[2/9] Calculando métricas de QC...")
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(adata.obs['n_genes_by_counts'], bins=60, color='steelblue', edgecolor='white')
axes[0].set_title('Genes por célula')
axes[1].hist(adata.obs['total_counts'], bins=60, color='salmon', edgecolor='white')
axes[1].set_title('UMIs por célula')
axes[2].hist(adata.obs['pct_counts_mt'], bins=60, color='mediumseagreen', edgecolor='white')
axes[2].axvline(x=20, color='red', linestyle='--')
axes[2].set_title('% Mitocondrial')
plt.tight_layout()
plt.savefig('qc_metrics.png', dpi=150, bbox_inches='tight')
plt.close()
print("  -> qc_metrics.png salvo")

 # --- 3. Filtragem ---
print("\n[3/9] Filtrando células e genes...")
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs['n_genes_by_counts'] < 5000, :]
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
print(f"  {adata.n_obs} células, {adata.n_vars} genes (após filtragem)")

# --- 4. Normalização ---
print("\n[4/9] Normalizando e aplicando log1p...")
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# --- 5. HVGs ---
print("\n[5/9] Selecionando genes altamente variáveis...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
n_hvg = adata.var['highly_variable'].sum()
print(f"  {n_hvg} HVGs selecionados")
adata = adata[:, adata.var['highly_variable']]

# --- 6. Regressão e escalonamento ---
print("\n[6/9] Regressando variáveis de confusão e escalando...")
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# --- 7. PCA ---
print("\n[7/9] Calculando PCA...")
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

# --- 8. Vizinhança e UMAP ---
print("\n[8/9] Construindo grafo e UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')

sc.pl.umap(adata, color='leiden_0.5', legend_loc='on data',
         title='Clusters Leiden (res=0.5)', save='_clusters.png')
print("  -> figures/umap_clusters.png salvo")

# --- 9. Marcadores ---
print("\n[9/9] Identificando genes marcadores...")
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', method='wilcoxon', use_raw=True)
markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_df.to_csv('markers_por_cluster.csv', index=False)
print("  -> markers_por_cluster.csv salvo")

# Verificar marcadores canônicos
marcadores = ['CD3E','CD4','CD8A','CD8B','IL7R','FOXP3',
              'CD19','MS4A1','CD79A','GNLY','NKG7',
              'CD14','LYZ','CD68','EPCAM','KRT18','KRT19']
presentes = [g for g in marcadores if g in adata.raw.var_names]
if presentes:
   sc.pl.umap(adata, color=presentes[:9], use_raw=True, ncols=3, save='_canonical_markers.png')
   print(f"  -> umap_canonical_markers.png salvo ({len(presentes)} marcadores encontrados)")

 # Salvar
adata.write('nsclc_scrna_analisado.h5ad')

print("\n" + "=" * 50)
print("  ANÁLISE CONCLUÍDA!")
print(f"  Células finais: {adata.n_obs}")
print(f"  Clusters: {adata.obs['leiden_0.5'].nunique()}")
print("  Arquivos gerados:")
print("    - qc_metrics.png")
print("    - figures/ (UMAPs e heatmaps)")
print("    - markers_por_cluster.csv")
print("    - nsclc_scrna_analisado.h5ad")
print("=" * 50)
