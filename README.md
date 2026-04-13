# BINF_6110_Assignment_4
scRNA-seq
## Introduction

Influenza A virus (IAV) is a negative-sense, single-stranded RNA virus of the family Orthomyxoviridae and a leading cause of acute respiratory illness, responsible for an estimated 290,000-650,000 deaths annually [1]. The nasal mucosa is the primary site of IAV entry and initial immune sensing. The murine nasal cavity comprises three anatomically and functionally distinct regions: the olfactory mucosa (OM), the respiratory mucosa (RM), and the lung-proximal nasal tissue (LNG), each containing unique epithelial and immune cell populations. Characterising how these compartments respond to infection across time is important for understanding mucosal antiviral defence and the establishment of long-lived respiratory immune memory.

IAV is detected by pattern recognition receptors including TLR3, TLR7, RIG-I and MDA5, which drive production of type I and III interferons that induce interferon-stimulated genes (ISGs) such as Isg15, Rsad2, and Ifit1, restricting viral replication [2]. NK cells and macrophages are recruited early, providing cytotoxic clearance and chemokine-mediated lymphocyte recruitment respectively [3]. Natural killer T (NKT) cells bridge innate and adaptive immunity through rapid cytokine production upon activation [4]. Following clearance, tissue-resident memory T cells (TRM) retained within mucosal niches mediate accelerated protection upon rechallenge, making their maintenance a key correlate of immunity.

Kazer et al. (2024) generated a comprehensive atlas of 156,572 murine nasal mucosal cells across the OM, RM, and LNG at five timepoints spanning primary IAV infection (naive, D02, D05, D08, D14) [5]. A key finding was the identification of KNIIFE cells (Krt13+ nasal immune-interacting floor epithelial cells), defined by co-expression of Krt13, Il1a, Cd274 (PD-L1), and Cxcl16. KNIIFE cells reside along the anterior nasal floor and expand following viral clearance at D14. Their expression of CXCL16, the sole ligand for CXCR6, implicates them as niche-providing cells for CXCR6+ tissue-resident CD8+ T cells, suggesting a role in long-term immune memory maintenance [5]. Co-expression of PD-L1 further implies a tolerogenic function limiting post-infectious immunopathology.

Single-cell RNA sequencing (scRNA-seq) is well-suited to characterising complex mucosal tissues during infection. Unlike bulk RNA-seq, which averages signals across heterogeneous populations, scRNA-seq captures transcriptome-wide profiles at single-cell resolution, enabling unbiased identification of rare populations and cell-type-specific transcriptional responses [6]. Analysis was performed using the Seurat v4 framework [7], which provides an integrated pipeline for quality control, normalisation, dimensionality reduction, clustering, and differential expression.

Harmony [8] was selected for batch correction over Seurat CCA [7], BBKNN [9], and Scanorama [10], as it operates on the PCA embedding to efficiently remove batch variance while preserving biological structure, and avoids the overcorrection artefacts associated with CCA. Fourteen principal components were retained based on elbow plot inflection [7]. For cell type annotation, SingleR [11] was applied against the ImmGen reference [12], which offers granular murine immune cell profiles and quantitative confidence scores, outperforming Azimuth (human-focused) and scType (limited tissue resolution) for this dataset. A manual KNIIFE annotation was performed using an additive module score of four canonical markers (Krt13, Il1a, Cd274, Cxcl16), following Kazer et al. [5]. Differential expression was performed using the Wilcoxon rank-sum test, as insufficient biological replicates per group blocked pseudobulk methods such as DESeq2 [13] or edgeR [14], with the acknowledged limitation that treating cells as independent observations may inflate false positive rates [15]. Pathway enrichment was assessed using GSEA [16] against the MSigDB Hallmark [17] and C7 immunological signature collections, and ORA against Gene Ontology Biological Process terms [18] for cluster marker gene sets.

This study applies the above pipeline to the Kazer et al. (2024) primary infection dataset with four objectives: (1) to perform quality-controlled preprocessing and Harmony-integrated clustering across five timepoints and three nasal regions; (2) to annotate cell populations using KNIIFE module scoring and SingleR; (3) to characterise the transcriptional dynamics of a neutrophil-enriched cluster (cluster 5) and an NKT cell-enriched cluster (cluster 6); and (4) to quantify the innate immune response of RM NK cells and macrophages at peak infection (D08) and contextualise findings within the mucosal antiviral immunity literature.

## Methods

### Dataset

The processed single-cell RNA sequencing dataset was obtained from Kazer et al. (2024) [5], accessed via GEO accession GSE266469 and the Broad Single Cell Portal (SCP2216). The dataset comprises murine nasal mucosal cells collected from three anatomical regions (olfactory mucosa [OM], respiratory mucosa [RM], and lung-proximal tissue [LNG]) across five timepoints during primary influenza A virus infection: naive, day 2 (D02), day 5 (D05), day 8 (D08), and day 14 (D14) post-infection. Data were provided as a pre-processed Seurat object. Sample identity was encoded in the orig.ident metadata field, from which timepoint and tissue region were parsed using regular expressions. All analyses were conducted in R v4.3 [19] using Seurat v4 [7]. A global random seed of 42 was set at the start of the analysis to ensure reproducibility.

### Quality Control

Quality control metrics were calculated per cell: the number of detected genes (nFeature_RNA), total UMI count (nCount_RNA), and percentage of reads mapping to mitochondrial genes (percent.mt), computed using PercentageFeatureSet() with the pattern ^mt-. Cells were retained if they satisfied all of the following criteria: 500 < nFeature_RNA < 2,500; 750 < nCount_RNA < 6,000; and percent.mt < 8%. Lower thresholds on gene and UMI counts excluded low-quality or empty droplets, upper thresholds excluded potential doublets, and the mitochondrial threshold excluded cells with compromised membrane integrity. These thresholds were applied using subset() and QC metrics were visualised with violin plots before and after filtering.

### Normalisation, Cell Cycle Scoring, and Feature Selection

Counts were normalised using log-normalisation (NormalizeData(), method = "LogNormalize", scale factor = 10,000), which divides each cell's counts by its total count, multiplies by the scale factor, and applies a natural log transformation. Cell cycle phase scores (S and G2M) were assigned using CellCycleScoring() with the canonical Seurat cell cycle gene lists, converted to title case for compatibility with mouse gene nomenclature using str_to_title(). The top 3,000 highly variable features were identified using variance-stabilising transformation (FindVariableFeatures(), selection.method = "vst"). Data were then scaled across all cells using ScaleData(), regressing out mitochondrial content (percent.mt) and cell cycle scores (S.Score, G2M.Score) to remove confounding sources of variation from the embedding.

### Dimensionality Reduction and Batch Correction

Principal component analysis (PCA) was performed on the scaled variable features using RunPCA(). The number of PCs to retain was determined by visual inspection of the elbow plot, with 14 PCs selected at the point of inflection in the standard deviation curve. Batch correction across samples was performed using Harmony [8] (RunHarmony(), group.by.vars = "orig.ident", dims.use = 1:14), which iteratively adjusts PC embeddings to minimise sample-driven variance while preserving biological structure. UMAP visualisation was computed on the Harmony-corrected embedding (RunUMAP(), reduction = "harmony", dims = 1:14, n.neighbors = 30). Batch correction quality was assessed by inspecting the UMAP coloured by sample identity.

### Graph-Based Clustering

A shared nearest-neighbour (SNN) graph was constructed on the Harmony embedding using FindNeighbors() (reduction = "harmony", dims = 1:14). Clusters were identified using the Louvain algorithm implemented in FindClusters() at resolution 0.3, a value selected to yield a biologically interpretable number of clusters without oversegmentation. Cluster identity was stored in the seurat_clusters metadata column. Canonical marker gene expression across clusters was visualised using dot plots to support cell type characterisation.

### Cell Type Annotation

Two annotation strategies were applied. First, a manual approach identified KNIIFE cells using AddModuleScore() with four canonical markers (Krt13, Il1a, Cd274, Cxcl16), following Kazer et al. [5]. Clusters with a mean module score >= 0.20 were designated KNIIFE; all others were labelled Other. This threshold was chosen based on the score distribution to selectively capture clusters with robust co-expression of all four markers.
Second, automated annotation was performed using SingleR [11] against the ImmGen reference dataset (celldex::ImmGenData()) [12], using main-level labels and Wilcoxon-based marker detection (de.method = "wilcox"). SingleR assigns each cell a label by computing Spearman correlations with reference profiles and applies a delta-score pruning step to flag low-confidence annotations. The resulting pruned labels (singler_pruned) were used for downstream compositional analysis. Annotation confidence was assessed using score heatmaps (plotScoreHeatmap()) and delta score distributions (plotDeltaDistribution()).

### Differential Expression Analysis

Differential expression (DE) analysis was performed using FindMarkers() in Seurat [7]. Positive cluster marker genes for clusters 5 and 6 were identified against all other clusters using the bimodal likelihood ratio test (test.use = "bimod") [20], with minimum cell fraction thresholds of min.pct = 0.35 and log2 fold-change threshold of 0.5. Timepoint-specific DE was performed using the Wilcoxon rank-sum test: cluster 5 cells were compared between D14 and naive, and cluster 6 cells between D05 and naive (min.pct = 0.15, logfc.threshold = 0.20). An additional Wilcoxon-based comparison was performed between D08 and naive within RM-resident NK cells and macrophages (identified by SingleR labels), using the same thresholds. A pseudobulk approach (DESeq2 [13]) was not feasible due to insufficient biological replicates per condition, which is acknowledged as a limitation. All p-values were adjusted using the Bonferroni correction as implemented in Seurat.

### Pathway Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) was performed using the GSEA() function from the clusterProfiler package [21] (v4.8). Ranked gene lists were constructed from the average log2 fold-change values from each DE comparison. Gene sets were retrieved using the msigdbr package [22] for Mus musculus, querying the MSigDB Hallmark collection (category "H") [17] and the C7 immunological signatures collection (category "C7") for the innate immune comparison. Significance threshold was set at adjusted p-value < 0.05. Over-representation analysis (ORA) using Gene Ontology Biological Process (GO-BP) terms [18] was performed with enrichGO() for upregulated marker genes from clusters 5 and 6 (avg_log2FC > 0.5, p_val_adj < 0.05), using all expressed genes in the dataset as the background universe. Gene symbols were mapped to Entrez IDs via bitr() against the org.Mm.eg.db annotation package. Multiple testing correction used the Benjamini-Hochberg method (p-value cutoff = 0.05, q-value cutoff = 0.2).

### Cell Composition Analysis

Cell type composition across timepoints and tissues was quantified using SingleR pruned labels. Fractional abundance of each cell type was calculated per timepoint within RM cells and per tissue at D08. Krt13 expression kinetics across infection were summarised as mean log-normalised expression per timepoint in RM cells, computed using FetchData(), serving as a proxy for KNIIFE cell activity over the course of infection. Composition was visualised using stacked bar plots with ggplot2 [23].

## Results

### Quality Control and Dataset Overview

After applying quality control filters (500 < nFeature_RNA < 2,500; 750 < nCount_RNA < 6,000; percent.mt < 8%), 97,747 cells were retained across 25,129 genes. Post-filtering, cells had a median of 1,568 detected genes (range 505-2,499), a median UMI count of 4,408 (range 3,247-5,999), and a median mitochondrial fraction of 2.75% (range 0-8.0%), indicating high-quality transcriptomic profiles (Figure 1). The five timepoints were represented by 23,927 naive, 24,197 D02, 11,736 D05, 16,094 D08, and 21,793 D14 cells.

<img width="975" height="348" alt="image" src="https://github.com/user-attachments/assets/b6650896-007d-4c0c-800b-cc20df0a2b0a" />

Figure 1: QC metrics after filtering. Violin plots showing the per-cell distributions of nFeature_RNA, nCount_RNA, and percent.mt across all five infection timepoints following application of quality control thresholds. Cells were retained if they had between 500 and 2,500 detected genes, between 750 and 6,000 UMIs, and a mitochondrial fraction below 8%. These filters removed low-quality, empty, and potentially doublet-containing droplets. Post-filtering distributions are consistent and well-centred across all timepoints, indicating uniform data quality.

### Dimensionality Reduction, Batch Correction, and Clustering

PCA was performed on the 3,000 most variable genes. Inspection of the elbow plot revealed an inflection at 14 principal components, which were retained for downstream analysis (Figure 2). Harmony integration across the 15 samples produced a well-mixed UMAP embedding with no visually apparent sample-specific clustering (Figure 3), confirming effective batch correction. Graph-based clustering at resolution 0.3 yielded 20 transcriptionally distinct clusters (Figure 4), with sizes ranging from 79 (cluster 19) to 18,446 (cluster 0) cells. UMAP visualisation by tissue region (Figure 5) and timepoint (Figure 6) revealed broad intermingling of regions and timepoints across clusters, consistent with shared cell type identities across anatomical sites and infection states

<img width="634" height="363" alt="image" src="https://github.com/user-attachments/assets/d1233bf9-4b59-4d79-9afe-f02688d285f5" />

Figure 2: Elbow plot of standard deviation explained by each principal component. The curve shows a clear inflection between PC 10 and PC 15, with the red dashed line marking the selected cutoff of 14 PCs. Retaining 14 PCs balances biological signal capture against the inclusion of uninformative dimensions in the downstream Harmony integration and neighbour graph construction.


<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/f72c213c-afd1-4256-a130-4eaca56b862d" />

Figure 3: UMAP coloured by sample identity (orig.ident), used to assess the effectiveness of Harmony batch correction. Each colour represents one of the 15 samples spanning three tissue regions and five timepoints. The well-mixed distribution of all samples across every cluster, with no visible sample-specific enrichment, confirms that Harmony successfully removed technical batch-driven variance while preserving the underlying biological cell type structure. Samples from all timepoints and tissues co-localise within the same cluster regions.

<img width="634" height="555" alt="image" src="https://github.com/user-attachments/assets/e27ba771-a8c8-4d43-91bc-b49eeb6d9b67" />

Figure 4: UMAP coloured by Seurat cluster identity. Twenty transcriptionally distinct clusters were identified using graph-based Louvain clustering at resolution 0.3 on the Harmony-corrected PCA embedding. Cluster sizes ranged from 79 cells (cluster 19) to 18,446 cells (cluster 0). Cluster labels are displayed at centroid positions. The spatial separation of clusters on the UMAP reflects underlying transcriptional distance, with closely positioned clusters sharing greater gene expression similarity.

<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/84ec9e8e-db45-44e5-af6b-bfa23ce194a7" />

Figure 5: UMAP coloured by anatomical tissue region. The three regions are broadly distributed across clusters rather than segregating into region-specific groupings, indicating that the dominant axis of transcriptional variation across this dataset is cell type identity rather than tissue of origin. Some clusters show moderate regional enrichment, reflecting tissue-specific cell type abundance differences.

<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/2de986aa-c6ee-4762-8060-beaf3db9977b" />

Figure 6: UMAP coloured by infection timepoint. Naive cells are shown in light grey; infected timepoints are coloured in a progression from red through orange, yellow, to green. Cells from all timepoints co-distribute across the major cluster regions, indicating that core cell type identities are preserved throughout infection. Subtle shifts in cell density within clusters between timepoints reflect infection-driven changes in cell type composition and transcriptional state.

### Cell Type Annotation

Canonical marker gene expression across all 20 clusters was visualised using a dot plot (Figure 7). Epithelial markers (Epcam, Pigr) were broadly expressed across multiple clusters, while Krt13 showed high specificity for clusters 1 and 18. Macrophage markers (C1qa, Csf1r) were confined to cluster 1, T/NKT cell markers (Cd3e, Cd8a) to cluster 6, and neutrophil markers (S100a8, Ly6c2) to cluster 5. Interferon-stimulated genes (Isg15, Ifitm3) were broadly expressed across multiple clusters, consistent with a systemic antiviral response.

<img width="975" height="488" alt="image" src="https://github.com/user-attachments/assets/9c98b994-a3df-4f76-bf4d-c7eb180100c0" />

Figure 7: Dot plot of canonical marker gene expression across all 20 Seurat clusters. Genes represent markers for epithelial cells (Krt13, Epcam, Pigr), NK cells (Nkg7, Gzmb, Prf1), macrophages (C1qa, Csf1r, Adgre1), monocytes (Ly6c2, S100a8), B cells (Ms4a1, Cd79a), T/NKT cells (Cd3e, Cd8a), and interferon-stimulated genes (Ifitm3, Isg15, Oas1a). Dot size indicates the percentage of cells within a cluster expressing each gene; colour intensity represents the average scaled expression level. The specificity of marker expression confirms the transcriptional identity of key clusters and validates the clustering resolution.

Manual KNIIFE annotation using an additive module score of four markers (Krt13, Il1a, Cd274, Cxcl16) identified clusters 1 (mean score 0.41) and 18 (mean score 1.69) as KNIIFE-positive, yielding 10,059 KNIIFE cells (10.3% of total) and 87,688 Other cells (Figure 8). KNIIFE cells were distributed across all three tissues (LNG: 4,211; OM: 2,245; RM: 3,603) and all five timepoints. Automated SingleR annotation against the ImmGen reference assigned 20 cell type labels, with epithelial cells (n = 21,662), fibroblasts (n = 13,497), and B cells (n = 12,722) constituting the largest populations (Figure 9). Only 542 cells (0.6%) received low-confidence pruned labels. Delta score distributions confirmed high annotation confidence across most cell types, with T cell lineages (T cells, Tgd, ILC) showing expectedly lower delta values due to transcriptional similarity (Figure 10). SingleR correctly assigned the KNIIFE-containing clusters as epithelial cells, consistent with the manual annotation.

<img width="731" height="506" alt="image" src="https://github.com/user-attachments/assets/2484eec6-3e23-4777-84ac-e958a56e4c0d" />

Figure 8: UMAP coloured by manual KNIIFE annotation. Cells are coloured orange or grey based on whether their parent cluster exceeded a mean module score threshold of 0.20 for the four-gene KNIIFE signature (Krt13, Il1a, Cd274, Cxcl16). Clusters 1 (score = 0.41) and 18 (score = 1.69) were classified as KNIIFE-positive, yielding 10,059 KNIIFE cells (10.3% of the dataset). KNIIFE cells occupy a discrete region of the UMAP, spatially separated from immune and stromal populations, consistent with their epithelial identity and transcriptional distinctiveness.

<img width="731" height="506" alt="image" src="https://github.com/user-attachments/assets/000be2d3-f6ea-4025-8761-d3d5e8d21a5e" />

Figure 9: UMAP coloured by SingleR automated cell type annotation against the ImmGen reference database. Twenty distinct cell type labels were assigned across 97,205 confidently annotated cells. The largest populations include epithelial cells, fibroblasts, B cells, dendritic cells, and macrophages. Grey cells represent low-confidence annotations pruned by the delta-score threshold. The spatial distribution of annotated populations is broadly consistent with the Seurat cluster structure, with immune populations co-localising in the lower-left region of the UMAP and epithelial and stromal populations occupying the right and upper regions.

<img width="780" height="467" alt="image" src="https://github.com/user-attachments/assets/a5bde7dc-c002-4eb1-8484-43a21ee48ee2" />

Figure 10: SingleR delta score distributions per annotated cell type. The delta score represents the difference between the top-scoring and second-scoring reference label for each cell, with higher values indicating greater annotation confidence and lower ambiguity. Most cell types display well-separated delta score distributions with high median values, confirming robust annotation. T cells, Tgd cells, and ILCs show comparatively lower delta scores, reflecting the transcriptional similarity between these lymphocyte subsets at the resolution of main-level ImmGen labels.

Feature plots confirmed spatial co-localisation of Nkg7 (NK cells), C1qa (macrophages), Ly6c2 (monocytes), and Ms4a1 (B cells) with their expected UMAP regions (Figure 12). Interferon-stimulated genes Ifitm3, Isg15, Oas1a, and Rsad2 were broadly expressed across multiple cell types (Figure 13). Epcam expression was broadly maintained in epithelial clusters at both naive and D14, confirming epithelial cluster stability across infection (Figure 15).


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f1676ee4-390c-4a5d-9d72-cc316b473619" />

######## delete ####### Figure 11: Violin plots of KNIIFE marker expression in KNIIFE vs Other cells.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/407e6a66-cdcd-437b-94dc-bc036dd1c392" />

Figure 12: Feature plots of immune cell marker genes projected onto the UMAP embedding. Nkg7, a cytotoxic granule protein, localises tightly to the upper-left cluster consistent with NK and cytotoxic T cell populations. C1qa, a complement component expressed by myeloid cells, shows highly specific expression in a single discrete cluster corresponding to macrophages. Ly6c2 marks a broader monocyte/granulocyte region spanning multiple clusters. Ms4a1 (CD20) localises to a compact cluster consistent with B cells. The spatial specificity of these markers validates the cell type assignments made by both manual inspection and SingleR annotation.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f4b06137-9521-4c8e-9851-3188de8d1ccf" />

Figure 13: Feature plots of four interferon-stimulated genes (Ifitm3, Isg15, Oas1a, Rsad2) projected onto the UMAP. In contrast to the cell-type-specific markers in Figure 12, ISGs are broadly expressed across multiple clusters spanning epithelial, immune, and stromal cell types, consistent with a systemic, cell-type-non-specific antiviral transcriptional response driven by paracrine interferon signalling during IAV infection. Expression is heterogeneous within clusters, reflecting variation in ISG induction levels between individual cells even within the same population.

<img width="780" height="520" alt="image" src="https://github.com/user-attachments/assets/d0044848-8cb6-44f1-8131-79c71e10f05d" />

####### delete ######## Figure 14: Violin plots of ISG expression (Ifitm3, Isg15, Oas1a, Rsad2) in KNIIFE vs Other cells. KNIIFE cells show elevated ISG expression, suggesting an interferon-responsive epithelial state.


<img width="828" height="414" alt="image" src="https://github.com/user-attachments/assets/06e141b0-226e-49da-85b7-24350428cfc8" />

Figure 15: Feature plots of Epcam expression split by timepoint, comparing naive (left) and day 14 post-infection (right) cells. Epcam, a pan-epithelial marker, is expressed broadly across the same epithelial cluster regions at both timepoints, confirming that the epithelial compartment is maintained structurally throughout the course of infection and into the resolution phase. The persistence of Epcam-positive clusters at D14 indicates that the mucosal epithelium recovers from IAV-induced damage and re-establishes normal epithelial architecture by the end of the primary infection time course.

### Cluster Marker Genes

FindAllMarkers identified positive marker genes for all 20 clusters (min.pct = 0.35, logFC > 0.5). The top five marker genes per cluster, visualised as a heatmap (Figure 16), demonstrated clear transcriptional separation between clusters. Cluster 5 was characterised by high expression of Retnlg (avg log2FC = 11.53), Cxcr2 (9.51), Mmp8 (8.56), Mmp9 (8.13), S100a9 (6.96), and S100a8 (6.75), consistent with a neutrophil identity [24]. In total, 342 of 345 tested genes were significant markers for cluster 5 (adjusted p < 0.05). Cluster 6 was defined by Cd3e (10.33), Cd3g (10.19), Trbc1 (10.16), Il2rb (9.99), Cxcr6 (9.90), and Gzmb (8.79), consistent with an NKT cell identity [4]. 432 of 434 tested genes were significant markers for cluster 6.


<img width="975" height="731" alt="image" src="https://github.com/user-attachments/assets/62982e6a-b08d-4f5b-95a3-3b70e39fecdb" />

Figure 16: Heatmap of the top 5 positive marker genes per Seurat cluster, identified by FindAllMarkers (bimod test, min.pct = 0.35, logFC threshold = 0.5). Scaled expression values are shown across all cells, grouped by cluster identity along the x-axis. Each row represents a gene and each column a cell. The diagonal pattern of enrichment confirms that each set of top marker genes is specific to its target cluster with minimal cross-cluster expression, validating the transcriptional distinctiveness of all 20 clusters. Key clusters are readily identifiable by their marker blocks, including clusters 5 (Retnlg, Mmp8/9), 6 (Cd3e, Trbc1), and 18 (Krt13).

### Temporal Dynamics of Cluster 5 (Neutrophils) and Cluster 6 (NKT Cells)

The constitutive identity markers of both clusters were stable across all five timepoints (Figure 17, 18), confirming that cluster membership reflects true cell type identity rather than infection-state artefacts. Timepoint-specific DE revealed divergent responses between the two populations. In cluster 5 neutrophils, comparison of D14 versus naive yielded 117 significantly upregulated and 152 downregulated genes. Upregulated genes at D14 included Isg15 (log2FC = 1.36), Tgm2 (1.32), and Tnfaip2 (1.25), suggesting a mild, sustained interferon-responsive state. Downregulated genes included chromatin-associated factors Hmgn1 (-1.33) and Hmgn3 (-1.31), indicating transcriptional reprogramming at the resolution phase of infection.

In contrast, cluster 6 NKT cells mounted a substantially stronger transcriptional response at D05 versus naive, with 1,472 upregulated and 208 downregulated genes. The response was dominated by ISGs: Oasl1 (log2FC = 9.27), Rsad2 (8.67), Ifit3b (6.78), Ifit1 (6.54), Cmpk2 (6.30), and Isg15 (4.44), reflecting a robust antiviral transcriptional programme at peak early infection. GSEA confirmed significant enrichment of HALLMARK_INTERFERON_ALPHA_RESPONSE (NES = 2.71) and HALLMARK_INTERFERON_GAMMA_RESPONSE (NES = 2.71) in cluster 6 at D05, alongside suppression of HALLMARK_OXIDATIVE_PHOSPHORYLATION (NES = -3.08) (Figure 19, 20). In cluster 5 at D14, GSEA identified suppression of oxidative phosphorylation (NES = -2.96) and MYC targets (NES = -2.33), with modest activation of HALLMARK_TNFA_SIGNALING_VIA_NFKB (NES = 1.92) and interferon pathways. ORA of cluster 5 marker genes against GO Biological Process terms revealed enrichment of neutrophil chemotaxis, granulocyte migration, and actin filament organisation (Figure 21), while cluster 6 markers were enriched for cytoplasmic translation, T cell differentiation, and lymphocyte proliferation (Figure 22), consistent with an actively proliferating cytotoxic lymphocyte population.

<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f2a2ac78-b3bf-4683-b16f-f243e63c7f82" />

Figure 17: Violin plots of the top four marker genes for cluster 5 (Retnlg, Cxcr2, Mmp8, Mmp9) across all five infection timepoints. Retnlg encodes resistin-like molecule gamma, a canonical neutrophil granule protein; Cxcr2 is the primary neutrophil chemokine receptor; and Mmp8 and Mmp9 are neutrophil-derived matrix metalloproteinases. All four markers are expressed at stable, consistently high levels across naive, D02, D05, D08, and D14, confirming that cluster 5 maintains a core neutrophil transcriptional identity throughout the course of infection and is not subject to dramatic infection-driven identity shifts.

<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/b85a4384-b501-42ea-be11-305bf80f2b65" />

Figure 18: Violin plots of the top four marker genes for cluster 6 (Cd3e, Cd3g, Trbc1, Il2rb) across all five infection timepoints. Cd3e and Cd3g are invariant components of the T cell receptor-CD3 signalling complex; Trbc1 encodes the T cell receptor beta constant chain; and Il2rb (CD122) is a component of the IL-2 and IL-15 receptor expressed on NKT and NK cells. All four markers are stably expressed across all timepoints, confirming that cluster 6 retains its NKT/cytotoxic T cell identity throughout infection. This stability contrasts with the strong transcriptional response observed in the timepoint DE analysis, indicating that infection-driven changes represent a superimposed activation state rather than a change in cell identity.


<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/ba8b917f-8968-4684-a13d-d5e665f59b02" />

Figure 19: GSEA MSigDB Hallmark dot plot for cluster 6 (NKT cells), comparing D05 versus naive. Each dot represents one significantly enriched gene set (adjusted p < 0.05); dot size indicates the number of leading-edge genes and colour indicates significance. Activated gene sets (left panel) are dominated by HALLMARK_INTERFERON_ALPHA_RESPONSE (NES = 2.71) and HALLMARK_INTERFERON_GAMMA_RESPONSE (NES = 2.71), consistent with the massive ISG upregulation observed at D05. HALLMARK_INFLAMMATORY_RESPONSE and HALLMARK_ALLOGRAFT_REJECTION are also activated. In the suppressed panel (right), HALLMARK_OXIDATIVE_PHOSPHORYLATION (NES = -3.08) and HALLMARK_ADIPOGENESIS are significantly downregulated, indicating metabolic reprogramming of NKT cells towards a glycolytic, effector-active state during peak infection.

<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/99e8ea75-e36f-4900-ba3b-a306ac653efd" />

Figure 20: GSEA MSigDB Hallmark dot plot for cluster 5 (neutrophils), comparing D14 versus naive. Activated gene sets include HALLMARK_TNFA_SIGNALING_VIA_NFKB (NES = 1.92), HALLMARK_INTERFERON_GAMMA_RESPONSE (NES = 1.98), HALLMARK_INTERFERON_ALPHA_RESPONSE (NES = 1.81), and HALLMARK_HYPOXIA (NES = 1.87), suggesting residual pro-inflammatory signalling and an adaptive response to the hypoxic tissue environment at the resolution phase of infection. Suppressed gene sets include HALLMARK_OXIDATIVE_PHOSPHORYLATION (NES = -2.96), HALLMARK_MYC_TARGETS_V1 (NES = -2.33), and HALLMARK_ADIPOGENESIS (NES = -2.65), indicating downregulation of anabolic and proliferative programmes in late-infection neutrophils relative to naive baseline.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f744f499-ab57-4761-b030-4231e9c06a8c" />

Figure 21: Over-representation analysis (ORA) of Gene Ontology Biological Process (GO-BP) terms for cluster 5 upregulated marker genes (avg_log2FC > 0.5, adjusted p < 0.05). Dot size represents the number of genes in the leading edge; colour intensity represents statistical significance. The top enriched terms are dominated by neutrophil and granulocyte biology: actin filament organisation, chemotaxis, leukocyte migration, granulocyte migration, neutrophil chemotaxis, and superoxide anion generation are all significantly enriched. The presence of neutrophil-specific GO terms at the highest gene ratios (0.10-0.13) strongly confirms the neutrophil identity of cluster 5 and is consistent with its marker gene profile including Retnlg, Cxcr2, Mmp8, Mmp9, S100a8, and S100a9.

<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/925ef457-2249-473b-8fb6-a836120738fc" />

Figure 22: ORA of GO-BP terms for cluster 6 upregulated marker genes. The most significantly enriched term is cytoplasmic translation (gene ratio ~0.20, adjusted p ~1e-65), reflecting the high ribosomal and translational activity of an actively effector-committed lymphocyte population. Additional enriched terms include ribonucleoprotein complex biogenesis, ribosome biogenesis, T cell differentiation, alpha-beta T cell activation, lymphocyte differentiation, lymphocyte proliferation, and leukocyte cell-cell adhesion. Together these terms are highly consistent with an NKT or cytotoxic T cell identity characterised by active protein synthesis, immune cell activation, and proliferative capacity.

### Innate Immune Response in RM NK Cells and Macrophages at D08

To characterise the innate immune response at peak infection, differential expression was performed between D08 and naive within RM-resident NK cells and macrophages identified by SingleR. This comparison yielded 1,226 significantly upregulated and 220 downregulated genes (adjusted p < 0.05). The most strongly upregulated genes included Cxcl9 (log2FC = 8.08), Gbp2 (5.16), Ifi205 (5.14), Rsad2 (4.49), Cd274 (4.21), and Isg15 (3.65). Downregulated genes included S100a5 (-3.73), Omp (-2.58), and Stoml3 (-2.50), reflecting suppression of baseline tissue-homeostatic programmes (Figure 23).
GSEA of the ranked gene list against the MSigDB Hallmark collection identified seven significantly enriched gene sets (Figure 24, 25). The two most strongly activated pathways were HALLMARK_INTERFERON_GAMMA_RESPONSE (NES = 2.57, adjusted p = 2.15 x 10-9) and HALLMARK_INTERFERON_ALPHA_RESPONSE (NES = 2.54, adjusted p = 2.15 x 10-9), reflecting the dominant role of interferon signalling in shaping innate immune cell states at peak infection. HALLMARK_ALLOGRAFT_REJECTION (NES = 1.96) and HALLMARK_INFLAMMATORY_RESPONSE (NES = 1.91) were also significantly activated, indicating broad immune effector and pro-inflammatory programmes. Conversely, HALLMARK_OXIDATIVE_PHOSPHORYLATION was suppressed (NES = -1.94), consistent with a metabolic shift away from mitochondrial respiration during acute inflammatory activation. GSEA against the C7 immunological signatures (Figure 26) yielded multiple significantly enriched gene sets related to T cell and mast cell activation states, with the highest NES values (>2.8) associated with adenosine receptor-mediated immune activation programmes.

<img width="731" height="569" alt="image" src="https://github.com/user-attachments/assets/2c732e91-48f7-4e11-bd88-f795b21a981d" />

Figure 23: Volcano plot of differential expression in RM-resident NK cells and macrophages comparing D08 versus naive. Each point represents one gene; the x-axis shows average log2 fold-change and the y-axis shows -log10 adjusted p-value. Significantly upregulated genes at D08 are shown in red and downregulated genes in blue; non-significant genes are grey. Dashed lines indicate the fold-change (±0.25) and significance (adjusted p = 0.05) thresholds. Key labelled genes include the top upregulated chemokines Cxcl9 and Cxcl10, the interferon-stimulated genes Isg15 and Ifit1, and the pro-inflammatory cytokine Tnf. The strong asymmetry towards upregulation (1,226 up vs 220 down) reflects the dominant transcriptional activation of the innate antiviral response at peak infection.

<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/c744af20-7cb3-4666-85d4-dc9ea5599a39" />

Figure 24: GSEA MSigDB Hallmark dot plot for RM NK cells and macrophages comparing D08 versus naive. Each dot represents one significantly enriched Hallmark gene set (adjusted p < 0.05). Dot size indicates the number of leading-edge genes contributing to enrichment; colour indicates adjusted p-value. The two most strongly activated pathways are HALLMARK_INTERFERON_ALPHA_RESPONSE (NES = 2.54) and HALLMARK_INTERFERON_GAMMA_RESPONSE (NES = 2.57), both at adjusted p = 2.15 x 10-9, reflecting the dominant role of type I and II interferon signalling in shaping innate immune cell states at D08. HALLMARK_ALLOGRAFT_REJECTION (NES = 1.96) and HALLMARK_INFLAMMATORY_RESPONSE (NES = 1.91) are also activated. HALLMARK_OXIDATIVE_PHOSPHORYLATION is the sole suppressed pathway (NES = -1.94), consistent with a metabolic shift towards glycolysis during acute inflammatory activation.

<img width="828" height="552" alt="image" src="https://github.com/user-attachments/assets/d91f3cfb-94a8-4c11-ab74-b6a57af8cb20" />

Figure 25: Ridge plot of GSEA enrichment distributions for MSigDB Hallmark gene sets in RM NK cells and macrophages, D08 versus naive. Each ridge represents the distribution of enrichment scores for the leading-edge genes of one gene set. IFN-gamma and IFN-alpha response distributions are shifted strongly to the right, confirming robust positive enrichment. Allograft rejection and inflammatory response show broader, bimodal distributions indicative of heterogeneous enrichment across the gene set. IL6-JAK-STAT3 and IL2-STAT5 signalling show modest positive shifts. Oxidative phosphorylation is shifted to the left, with a tight, left-centred distribution consistent with coordinated downregulation of mitochondrial respiratory chain gene expression during the antiviral activation state.

<img width="975" height="488" alt="image" src="https://github.com/user-attachments/assets/647ca5e5-1cc2-4784-87d6-2ff9a97cdab8" />

Figure 26: GSEA C7 immunological signatures dot plot for RM NK cells and macrophages comparing D08 versus naive. The C7 collection contains gene sets derived from published immunological experiments, enabling benchmarking of the observed transcriptional response against established immune activation states. The most significantly enriched activated signatures (NES > 2.5, adjusted p < 5.64 x 10-9) include multiple adenosine A3 receptor-related T cell and mast cell activation gene sets, as well as CD4 vs CD8 T cell comparisons and Th1 vs Th17 polarisation signatures. The suppressed panel includes naive vs activated CD8 T cell signatures and several Treg-associated gene sets, suggesting that the D08 innate cell transcriptional state is most similar to that of activated effector immune cells rather than resting or regulatory populations.

### Cell Type Composition Across Timepoints and Tissues

Analysis of RM cell type composition across primary infection (Figure 27) revealed that epithelial cells and fibroblasts constituted the dominant fractions at naive baseline (29.1% and 21.3%, respectively). During peak infection (D02-D08), the relative abundance of immune populations including monocytes, macrophages, NKT cells, and neutrophils increased, reflecting immune cell infiltration and expansion. Epithelial and stromal fractions partially recovered by D14. Cross-tissue composition at D08 (Figure 28) showed broadly similar patterns across LNG, OM, and RM, with epithelial cells consistently the most abundant fraction in all three regions (LNG: 30.8%; RM and OM comparable). Macrophage abundance was proportionally higher in LNG (12.9%) relative to the other regions. Krt13 mean expression in RM cells peaked modestly at D05 (mean log-normalised expression = 0.034) before declining at D08 and D14 (Figure 29), suggesting transient induction of KNIIFE-associated transcriptional activity at early-to-mid infection rather than a sustained D14 expansion, which may reflect limitations of the cluster-level annotation approach used.

<img width="828" height="497" alt="image" src="https://github.com/user-attachments/assets/d3de36ab-cbfb-4a85-a93b-c9c6c8bd0306" />

Figure 27: Stacked bar chart of cell type composition in RM cells across all five primary infection timepoints, based on SingleR pruned label annotations. Each bar represents the fractional abundance of each annotated cell type within the RM at one timepoint. At naive baseline, epithelial cells (29.1%) and fibroblasts (21.3%) dominate, reflecting the structural composition of the uninfected respiratory mucosa. During active infection (D02-D08), the relative proportions of innate immune populations including monocytes, macrophages, neutrophils, and NKT cells increase, consistent with immune infiltration driven by viral replication. By D14, epithelial and stromal fractions begin to recover towards baseline proportions, indicating partial mucosal restitution following viral clearance.


<img width="731" height="548" alt="image" src="https://github.com/user-attachments/assets/8c191c7c-4def-4bcc-b5c1-449522dcbca7" />

Figure 28: Stacked bar chart of cell type composition across the three nasal tissue regions (LNG, OM, RM) at D08, the peak infection timepoint. Epithelial cells are the most abundant fraction in all three regions (LNG: 30.8%, with OM and RM comparable). Fibroblasts are the second most abundant in all regions. Macrophages are proportionally more abundant in LNG (12.9%) relative to OM and RM, potentially reflecting greater viral replication and tissue damage in the lung-proximal compartment. Monocytes are most abundant in RM (data shown), consistent with preferential monocyte recruitment to the respiratory mucosa. The broadly similar compositional profiles across tissues at D08 suggest that the peak innate immune response is a pan-nasal mucosal phenomenon rather than being restricted to a single anatomical compartment.

<img width="731" height="406" alt="image" src="https://github.com/user-attachments/assets/86abb63a-fafc-4cd3-8182-028f7c06dc09" />

Figure 29: Line plot of mean log-normalised Krt13 expression in RM cells across the five primary infection timepoints, used as a proxy for KNIIFE cell transcriptional activity. Expression is low at naive baseline (mean = 0.016), rises modestly at D02 (0.017) and peaks at D05 (0.034), before declining at D08 (0.013) and D14 (0.013). The transient peak at D05 rather than the sustained D14 expansion reported by Kazer et al. (2024) likely reflects the limitations of the cluster-level annotation used here, which captures all Krt13-expressing epithelial cells rather than isolating true KNIIFE cells at single-cell resolution. Nevertheless, the D05 peak is consistent with early induction of KNIIFE-associated transcriptional programmes during active viral replication.


## Discussion




## References







