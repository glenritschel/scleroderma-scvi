# SSc single-cell drug repurposing — final bundle

_Generated: 2025-09-18T15:38:44_

## What this is
- Prioritized compounds from LINCS ‘reversal’ of SSc skin scRNA-seq cluster signatures.
- Ranking combines overall reversal, fibrolineage selectivity, LINCS support breadth, and simple target presence checks.

## Top 10 (composite priority)

| base_compound   | moa                       | status         |   priority_final |   max_rev |   neglog10_best_p |
|:----------------|:--------------------------|:---------------|-----------------:|----------:|------------------:|
| PD-0325901      | MEK inhibitor             | Clinical-stage |        3.82521   |   28.0899 |           28.0899 |
| gefitinib       | EGFR inhibitor            | FDA-approved   |        2.53258   |   23.7427 |           23.7427 |
| WZ-3105         | Unknown/other             | Unknown        |        0.766407  |   15.3522 |           15.3522 |
| PD-184352       | MEK inhibitor             | Clinical-stage |        0.553922  |   22.2944 |           22.2944 |
| BI-2536         | PLK1 inhibitor            | Clinical-stage |        0.384485  |   21.4354 |           21.4354 |
| I-BET151        | BET bromodomain inhibitor | Unknown        |        0.275805  |   18.2887 |           18.2887 |
| JNK-9L          | JNK inhibitor             | Unknown        |        0.0970953 |   15.1053 |           15.1053 |
| I-BET           | BET bromodomain inhibitor | Clinical-stage |       -0.16189   |   18.8555 |           18.8555 |
| sirolimus       | mTOR inhibitor            | FDA-approved   |       -0.174531  |   18.4636 |           18.4636 |
| pelitinib       | EGFR inhibitor            | Unknown        |       -0.437989  |   17.0938 |           17.0938 |

## Files included
- final_shortlist_ranked.csv
- shortlist_shareable.csv
- dose_time_consistency_summary.csv
- lincs_reversal_top15_by_cluster.csv
- lincs_reversal_crosscluster.csv
- lincs_reversal_aggregate.csv
- candidate_slate_all.csv
- candidate_slate_fibrolineage.csv
- shortlist_fibro_clinical_view.csv
- shortlist_fibro_clinical_view_strict.csv
- shortlist_fibro_with_target_presence.csv
- shortlist_fibro_with_target_presence_raw.csv
- fig_final_shortlist_bar.png
- onepager_BI-2536.png
- onepager_I-BET.png
- onepager_I-BET151.png
- onepager_JNK-9L.png
- onepager_PD-0325901.png
- onepager_PD-184352.png
- onepager_WZ-3105.png
- onepager_gefitinib.png
- onepager_pelitinib.png
- onepager_sirolimus.png

## Notes & caveats
- ‘Target missing’ often reflects HVG/ID issues (e.g., Ensembl vs symbols) rather than true absence.
- Reversal scores come from LINCS cell lines; treat as hypotheses, not clinical guidance.
