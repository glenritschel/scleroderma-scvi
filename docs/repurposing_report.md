A cell-type–resolved transcriptomic reversal screen identifies candidate mechanisms for dermal fibrosis in systemic sclerosis
Abstract

We assembled a single-cell atlas of SSc skin and used LINCS L1000 transcriptomic reversal to nominate compounds that may counteract cell-type–specific disease programs, with emphasis on fibroblast and myofibroblast states. After QC, clustering, and curated annotation, we computed per-cluster differential expression and queried Enrichr/LINCS libraries. Aggregating reversal evidence across clusters yielded a candidate slate that was triaged by fibroblast selectivity, statistical support, and (when measurable) target presence in our data. The top mechanistic families emerging from this analysis include MEK/ERK pathway inhibition, EGFR axis inhibition, and BET bromodomain blockade, with broad cytotoxics (e.g., PLK1 or topo-II inhibitors) deprioritized for translational narratives. These findings require external validation and do not constitute treatment recommendations.

1. Data & preprocessing (summary)

Input: SSc skin scRNA-seq (GEO).

Workflow: scVI integration → neighborhood graph → Leiden clustering → UMAP.

QC: mitochondrial and ribosomal content re-computed; low-quality clusters flagged.

Annotation: panel-based scoring + DE sanity checks; manual curation produced the final set of labels (e.g., Fibroblasts, Myofibroblasts, Keratinocytes (basal/suprabasal), Endothelial, Pericyte/SMC, Lymphatic endothelium, T/NK, B/Plasma, Myeloid, DC, Mast).

Key figures/tables

UMAP with curated labels: results/figures/umap_celltype_provisional.png (or the latest from your notebook).

Cluster sizes: results/tables/cluster_sizes.csv.

DE table (Wilcoxon): results/tables/rank_genes_groups_leiden_wilcoxon.csv.

2. Reversal screen (methods)

We formed “disease signatures” per Leiden cluster from the DE table and queried Enrichr/LINCS:

Libraries: GO_Biological_Process_2023, Reactome_2022, KEGG_2021_Human, and LINCS_L1000_Chem_Pert_(up|down).

For each cluster, we extracted up/down gene lists (adj. p < 0.05 when available).

Reversal score: we used Enrichr output to quantify how strongly a compound’s perturbational signature opposes the cluster signature (captured in rev_score_* columns).

We then summarized at two levels:

Cluster level: Top 15 reversing compounds per cluster
→ results/drug_repurposing/lincs_reversal_top15_by_cluster.csv

Cross-cluster aggregation: Sum/max of reversal across clusters and number of distinct clusters hit
→ results/drug_repurposing/lincs_reversal_crosscluster.csv and …/lincs_reversal_aggregate.csv

3. Candidate ranking & triage

We assembled a candidate slate (candidate_slate_*.csv) and computed a composite priority balancing:

Total reversal across clusters (primary signal).

Fibrolineage selectivity (preference for fibroblast/myofibroblast clusters).

Statistical support (−log10 of the best p in the evidence).

Dose/time consistency where LINCS profiles offered multiple conditions
→ results/drug_repurposing/dose_time_consistency_summary.csv

Target presence/absence in our matrix (checked against raw/HVG; used for annotation/penalty only)
→ results/drug_repurposing/shortlist_fibro_with_target_presence_raw.csv

We also ran a weight sensitivity analysis to show that leading families are robust to the fibro-selectivity weight (figure saved earlier in your run).

Final ranked shortlist:
results/drug_repurposing/final_shortlist_ranked.csv
(Bar plot: results/drug_repurposing/final_shortlist_bar_top12.png)

Per-compound one-pagers:
results/drug_repurposing/final_bundle/onepager_*.png

4. Results (high-level)

Mechanistic families rising to the top (non-cytotoxic focus):

MEK/ERK axis: e.g., PD-0325901, PD-184352 (robust reversal; fibrolineage-relevant programs down).

EGFR axis: e.g., gefitinib (selective reversal in keratinocyte-rich/fibro-interactive contexts; merits cautionary interpretation as target detection in HVG/raw can be limited).

BET bromodomain: e.g., I-BET series (consistent, multi-cluster reversal; plausible anti-fibrotic transcriptional de-activation).

Cytotoxic “sweepers” (e.g., PLK1, topo-II) often score highly but are deprioritized for translational narratives due to lack of specificity.

Dose/time: where multiple LINCS conditions existed, top families showed coherent reversal across doses/timepoints (see dose_time_consistency_summary.csv and one-pagers).

Targets in data: direct target transcripts for some MOAs were not reliably detected in HVGs; we therefore report target presence checks using the raw base matrix and annotate “missing due to HVG/ID issues” where appropriate.

Illustrative artifacts

Heatmap of reversal (compound × cell type): see your earlier heatmap or browse lincs_reversal_top15_by_cluster.csv.

One-pager examples: onepager_PD-0325901.png, onepager_gefitinib.png, onepager_I-BET151.png.

5. Limitations

Model mismatch: LINCS signatures derive from cancer/immortalized lines that may not reflect SSc dermal fibroblasts or skin microenvironments.

Target detectability: HVG subsetting and gene ID formats (Ensembl vs symbols) can hide targets; we mitigated with raw checks but gaps remain.

No patient-level clinical stratification in this pass; results are atlas-level.

Reversal ≠ efficacy/safety: Transcriptomic opposition is hypothesis-generating only.

6. Next steps (Phase III)

External validation: rerun pipeline on an independent SSc skin dataset (single-cell or bulk fibroblast). Compare top MOAs and compound families.

Mechanistic sanity panels: quantify decrease of ECM/myofibroblast modules versus predicted reversal per cluster; include plots in the bundle.

Permutation null / empirical FDR: randomize gene labels and re-run 50–100× to estimate false discovery rates for the ranking.

In-vitro short list (research only): prioritize MEK/ERK, BET, and EGFR/mTOR families for dermal fibroblast assays (TGFβ-induced ECM; myofibroblast markers).

Heterogeneity: if sample labels exist, compute pseudo-bulk fibroblast signatures per sample and check whether leading MOAs generalize across samples.

Documentation: keep results/drug_repurposing/index.html and the final_bundle zip up to date; capture environment snapshots in results/metadata/.

7. Reproducibility

Artifacts: results/drug_repurposing/index.html links all CSVs/PNGs; final_bundle/ contains a shareable subset (shortlist CSV + one-pagers).

Environment: results/metadata/pip_freeze.txt and session.txt (Python/OS/git).

Key inputs: data/processed/ssc_skin_scvi_annot_curated.h5ad and results/tables/rank_genes_groups_leiden_wilcoxon.csv.

8. Ethical & clinical note

This work is a computational research screen. It does not assess dosing, interactions, or safety, and it should not be used to guide care. Any therapeutic investigation belongs in controlled preclinical models and clinical studies with physician oversight.