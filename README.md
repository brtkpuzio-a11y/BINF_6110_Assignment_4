# BINF_6110_Assignment_4
scRNA-seq
## Introduction

Influenza A virus (IAV) is a negative-sense, single-stranded RNA virus of the family Orthomyxoviridae and a leading cause of acute respiratory illness, responsible for an estimated 290,000-650,000 deaths annually [1]. The nasal mucosa is the primary site of IAV entry and initial immune sensing. The murine nasal cavity comprises three anatomically and functionally distinct regions: the olfactory mucosa (OM), the respiratory mucosa (RM), and the lung-proximal nasal tissue (LNG), each harbouring unique epithelial and immune cell populations. Characterising how these compartments respond to infection across time is essential for understanding mucosal antiviral defence and the establishment of long-lived respiratory immune memory.
IAV is detected by pattern recognition receptors including TLR3, TLR7, RIG-I and MDA5, which drive production of type I and III interferons that induce interferon-stimulated genes (ISGs) such as Isg15, Rsad2, and Ifit1, restricting viral replication [2]. NK cells and macrophages are recruited early, providing cytotoxic clearance and chemokine-mediated lymphocyte recruitment respectively [3]. Natural killer T (NKT) cells bridge innate and adaptive immunity through rapid cytokine production upon activation [4]. Following clearance, tissue-resident memory T cells (TRM) retained within mucosal niches mediate accelerated protection upon rechallenge, making their maintenance a key correlate of heterosubtypic immunity.
Kazer et al. (2024) generated a comprehensive atlas of 156,572 murine nasal mucosal cells across the OM, RM, and LNG at five timepoints spanning primary IAV infection (naive, D02, D05, D08, D14) [5]. A key finding was the identification of KNIIFE cells (Krt13+ nasal immune-interacting floor epithelial cells), defined by co-expression of Krt13, Il1a, Cd274 (PD-L1), and Cxcl16. KNIIFE cells reside along the anterior nasal floor and expand following viral clearance at D14. Their expression of CXCL16, the sole ligand for CXCR6, implicates them as niche-providing cells for CXCR6+ tissue-resident CD8+ T cells, suggesting a role in long-term immune memory maintenance [5]. Co-expression of PD-L1 further implies a tolerogenic function limiting post-infectious immunopathology.
Single-cell RNA sequencing (scRNA-seq) is well-suited to characterising complex mucosal tissues during infection. Unlike bulk RNA-seq, which averages signals across heterogeneous populations, scRNA-seq captures transcriptome-wide profiles at single-cell resolution, enabling unbiased identification of rare populations and cell-type-specific transcriptional responses [6]. Analysis was performed using the Seurat v4 framework [7], which provides an integrated pipeline for quality control, normalisation, dimensionality reduction, clustering, and differential expression.
Harmony [8] was selected for batch correction over Seurat CCA [7], BBKNN [9], and Scanorama [10], as it operates on the PCA embedding to efficiently remove batch variance while preserving biological structure, and avoids the overcorrection artefacts associated with CCA. Fourteen principal components were retained based on elbow plot inflection [7]. For cell type annotation, SingleR [11] was applied against the ImmGen reference [12], which offers granular murine immune cell profiles and quantitative confidence scores, outperforming Azimuth (human-focused) and scType (limited tissue resolution) for this dataset. A complementary manual KNIIFE annotation was performed using an additive module score of four canonical markers (Krt13, Il1a, Cd274, Cxcl16), following Kazer et al. [5]. Differential expression was performed using the Wilcoxon rank-sum test, as insufficient biological replicates per group precluded pseudobulk methods such as DESeq2 [13] or edgeR [14], with the acknowledged limitation that treating cells as independent observations may inflate false positive rates [15]. Pathway enrichment was assessed using GSEA [16] against the MSigDB Hallmark [17] and C7 immunological signature collections, and ORA against Gene Ontology Biological Process terms [18] for cluster marker gene sets.
This study applies the above pipeline to the Kazer et al. (2024) primary infection dataset with four objectives: (1) to perform quality-controlled preprocessing and Harmony-integrated clustering across five timepoints and three nasal regions; (2) to annotate cell populations using KNIIFE module scoring and SingleR; (3) to characterise the transcriptional dynamics of a neutrophil-enriched cluster (cluster 5) and an NKT cell-enriched cluster (cluster 6); and (4) to quantify the innate immune response of RM NK cells and macrophages at peak infection (D08) and contextualise findings within the mucosal antiviral immunity literature.

## Methods
## Results
<img width="975" height="348" alt="image" src="https://github.com/user-attachments/assets/b6650896-007d-4c0c-800b-cc20df0a2b0a" />

Figure 1: QC metrics after filtering. Violin plots showing distributions of nFeature_RNA, nCount_RNA, and percent.mt across all five timepoints following quality control filtering.


<img width="634" height="363" alt="image" src="https://github.com/user-attachments/assets/d1233bf9-4b59-4d79-9afe-f02688d285f5" />

Figure 2: Elbow plot of PCA standard deviation. The red dashed line indicates the selected cutoff of 14 principal components.


<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/f72c213c-afd1-4256-a130-4eaca56b862d" />

Figure 3: UMAP coloured by sample identity. Well-mixed distribution of all 15 samples across clusters confirms effective Harmony batch correction.


<img width="634" height="555" alt="image" src="https://github.com/user-attachments/assets/e27ba771-a8c8-4d43-91bc-b49eeb6d9b67" />

Figure 4: UMAP coloured by Seurat cluster identity. Twenty clusters were identified at resolution 0.3. Cluster labels are shown.


<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/84ec9e8e-db45-44e5-af6b-bfa23ce194a7" />

Figure 5: UMAP coloured by tissue region. The three regions are broadly distributed across clusters.


<img width="634" height="477" alt="image" src="https://github.com/user-attachments/assets/2de986aa-c6ee-4762-8060-beaf3db9977b" />

Figure 6: UMAP coloured by timepoint. Naive cells are shown in grey; infected timepoints (D02-D14) are shown in progressively warmer colours.


<img width="975" height="488" alt="image" src="https://github.com/user-attachments/assets/9c98b994-a3df-4f76-bf4d-c7eb180100c0" />

Figure 7: Dot plot of canonical marker gene expression across all 20 Seurat clusters. Dot size represents the percentage of cells expressing each gene; colour intensity represents average scaled expression.


<img width="731" height="506" alt="image" src="https://github.com/user-attachments/assets/2484eec6-3e23-4777-84ac-e958a56e4c0d" />

Figure 8: UMAP coloured by manual KNIIFE annotation. KNIIFE cells (orange) correspond to clusters 1 and 18, which crossed the module score threshold of 0.20.


<img width="731" height="506" alt="image" src="https://github.com/user-attachments/assets/000be2d3-f6ea-4025-8761-d3d5e8d21a5e" />

Figure 9: UMAP coloured by SingleR automated annotation against the ImmGen reference. Twenty cell types were assigned. Grey cells indicate low-confidence pruned labels (n = 542).


<img width="780" height="467" alt="image" src="https://github.com/user-attachments/assets/a5bde7dc-c002-4eb1-8484-43a21ee48ee2" />

Figure 10: SingleR delta score distributions per cell type. Higher delta values indicate greater annotation confidence.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f1676ee4-390c-4a5d-9d72-cc316b473619" />

Figure 11: Violin plots of KNIIFE marker expression in KNIIFE vs Other cells.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/407e6a66-cdcd-437b-94dc-bc036dd1c392" />

Figure 12: Feature plots of immune cell markers on UMAP. Nkg7 (NK cells), C1qa (macrophages), Ly6c2 (monocytes), and Ms4a1 (B cells) localise to their expected cluster regions.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f4b06137-9521-4c8e-9851-3188de8d1ccf" />

Figure 13: Feature plots of interferon-stimulated genes (Ifitm3, Isg15, Oas1a, Rsad2) on UMAP. Broad expression across multiple clusters reflects a systemic antiviral transcriptional response.


<img width="780" height="520" alt="image" src="https://github.com/user-attachments/assets/d0044848-8cb6-44f1-8131-79c71e10f05d" />

Figure 14: Violin plots of ISG expression (Ifitm3, Isg15, Oas1a, Rsad2) in KNIIFE vs Other cells. KNIIFE cells show elevated ISG expression, suggesting an interferon-responsive epithelial state.


<img width="828" height="414" alt="image" src="https://github.com/user-attachments/assets/06e141b0-226e-49da-85b7-24350428cfc8" />

Figure 15: Feature plots of Epcam expression split by timepoint (Naive vs D14). Epithelial cluster identity is maintained across infection.


<img width="975" height="731" alt="image" src="https://github.com/user-attachments/assets/62982e6a-b08d-4f5b-95a3-3b70e39fecdb" />

Figure 16: Heatmap of the top 5 marker genes per Seurat cluster. Scaled expression is shown across all cells, grouped by cluster identity. Rows represent genes, columns represent cells.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f2a2ac78-b3bf-4683-b16f-f243e63c7f82" />

Figure 17: Violin plots of cluster 5 top marker genes (Retnlg, Cxcr2, Mmp8, Mmp9) across infection timepoints. Expression is stable, confirming consistent neutrophil identity throughout infection.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/b85a4384-b501-42ea-be11-305bf80f2b65" />

Figure 18: Violin plots of cluster 6 top marker genes (Cd3e, Cd3g, Trbc1, Il2rb) across infection timepoints. Expression is stable, confirming consistent NKT cell identity throughout infection.


<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/ba8b917f-8968-4684-a13d-d5e665f59b02" />

Figure 19: GSEA Hallmark dot plot for cluster 6 (NKT cells), D05 vs naive. IFN-alpha and IFN-gamma responses are significantly activated; adipogenesis is suppressed.


<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/99e8ea75-e36f-4900-ba3b-a306ac653efd" />

Figure 20: GSEA Hallmark dot plot for cluster 5 (neutrophils), D14 vs naive. Oxidative phosphorylation and MYC targets are suppressed; TNFa/NF-kB and interferon signalling are modestly activated.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/f744f499-ab57-4761-b030-4231e9c06a8c" />

Figure 21: ORA of GO Biological Process terms for cluster 5 upregulated marker genes. Enrichment of neutrophil chemotaxis, granulocyte migration, and actin organisation confirms neutrophil identity.


<img width="731" height="584" alt="image" src="https://github.com/user-attachments/assets/925ef457-2249-473b-8fb6-a836120738fc" />

Figure 22: ORA of GO Biological Process terms for cluster 6 upregulated marker genes. Top terms include cytoplasmic translation, T cell differentiation, and lymphocyte proliferation, consistent with an NKT cell identity.


<img width="731" height="569" alt="image" src="https://github.com/user-attachments/assets/2c732e91-48f7-4e11-bd88-f795b21a981d" />

Figure 23: Volcano plot of differential expression in RM NK cells and macrophages, D08 vs naive. Significantly upregulated genes (red) and downregulated genes (blue) are shown; key cytotoxic and myeloid activation genes are labelled.


<img width="878" height="502" alt="image" src="https://github.com/user-attachments/assets/c744af20-7cb3-4666-85d4-dc9ea5599a39" />

Figure 24: GSEA MSigDB Hallmark dot plot for RM NK cells and macrophages, D08 vs naive. IFN-gamma and IFN-alpha responses dominate the activated side; oxidative phosphorylation is suppressed.


<img width="828" height="552" alt="image" src="https://github.com/user-attachments/assets/d91f3cfb-94a8-4c11-ab74-b6a57af8cb20" />

Figure 25: Ridge plot of GSEA enrichment distributions for MSigDB Hallmark gene sets. IFN-gamma and IFN-alpha show strong positive enrichment; oxidative phosphorylation shows a left-shifted distribution consistent with suppression.


<img width="975" height="488" alt="image" src="https://github.com/user-attachments/assets/647ca5e5-1cc2-4784-87d6-2ff9a97cdab8" />

Figure 26: GSEA C7 immunological signatures for RM NK cells and macrophages, D08 vs naive. Multiple T cell and mast cell activation signatures are significantly enriched.


<img width="828" height="497" alt="image" src="https://github.com/user-attachments/assets/d3de36ab-cbfb-4a85-a93b-c9c6c8bd0306" />

Figure 27: Stacked bar chart of cell type composition in RM across primary infection timepoints. Each bar represents the fractional abundance of SingleR-annotated cell types at each timepoint.


<img width="731" height="548" alt="image" src="https://github.com/user-attachments/assets/8c191c7c-4def-4bcc-b5c1-449522dcbca7" />

Figure 28: Cell type composition by tissue at D08. Stacked bars show fractional abundance of cell types in LNG, OM, and RM at peak infection.


<img width="731" height="406" alt="image" src="https://github.com/user-attachments/assets/86abb63a-fafc-4cd3-8182-028f7c06dc09" />

Figure 29: Mean Krt13 log-normalised expression in RM cells across timepoints, used as a proxy for KNIIFE cell activity. Expression peaks transiently at D05 before declining at later timepoints.





















## Discussion
## References

