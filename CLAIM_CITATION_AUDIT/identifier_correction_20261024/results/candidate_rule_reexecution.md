# Phase 4 renal-origin candidate freeze audit

Frozen: 2026-08-09

- DEG input: `CLAIM_CITATION_AUDIT/identifier_correction_20261024/results/corrected_DEG.csv`; SHA-256 `fdcc9ee72c89022f51ba6458f66c40a3a5e5bdaf59dc8c207ff77bb82bcf0a46`.
- Orthology input: `CLAIM_CITATION_AUDIT/identifier_correction_20261024/ortholog_rat_human_1to1_corrected.csv`; SHA-256 `07e2c86dc69ae6cf75a996793de798f79afaaf8bf1ce766057d17d29fe20f99c`.
- Renal DEGs passing the frozen rule: 1734.
- Valid one-to-one mapping rows: 1285.
- Unique frozen human candidates: 1285.
- Rat source symbols represented by multiple one-to-one mapping rows: 6; retained with an explicit row-level flag because each human target remains unique.
- AAD labels or AAD expression values used in candidate construction: **no**.
- Candidate role: primary renal-origin feature universe for downstream fold-internal modeling; not a causal gene set.

## Fail-closed checks

- Duplicate renal DEG symbols: rejected.
- One human symbol reached from multiple renal source genes: rejected.
- Non-one-to-one or confidence<1 mapping rows: excluded.
- Threshold drift from adj.P.Val<0.05 and |logFC|>0.5: rejected by tests.

