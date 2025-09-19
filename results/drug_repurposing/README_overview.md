# Drug-repurposing overview

*Updated:* 2025-09-18

This figure summarizes our LINCS reversal analysis for scleroderma skin single-cell states.

- **A.** Composite priority of candidate compounds (higher is better).
- **B.** Heatmap of max LINCS reversal per cell type for the top compounds.
- **C.** Selectivity of fibro/myofibro vs. other cell types plotted against the breadth of affected cell types.

![overview](fig_overview.png)

**Headline findings**
- MEK/ERK pathway inhibition (e.g., **PD-0325901**, **PD-184352**) ranks consistently high.
- **BET bromodomain** and **PLK1** inhibition emerge as additional mechanisms.
- Rankings are robust to selectivity weighting; dose/time metadata are limited for some compounds.

*Caveats:* LINCS cell contexts differ from primary skin cell states; some drug targets are absent from the HVG matrix (checked against `.raw` where possible).

<!-- PROVENANCE:BEGIN -->
## Provenance & reproducibility

- **Git:** tag `v0.1.5`, commit `03fdc2f`
- **AnnData used:** `n/a`
- **DE table:** `results/tables/rank_genes_groups_leiden_wilcoxon.csv`
- **LINCS libraries:** GO_Biological_Process_2023, Reactome_2022, KEGG_2021_Human, LINCS_L1000_Chem_Pert_up, LINCS_L1000_Chem_Pert_down
- **Environment:** Python 3.10.18; scanpy 1.11.4; scvi-tools 1.3.3; pandas 2.3.2; gseapy 1.1.10
- **Generated:** 2025-09-18 22:00:10

_Reproduce:_ run notebooks `06_results` → `07_state_signatures_and_drugs` → `08_validation_and_robustness`.
<!-- PROVENANCE:END -->
