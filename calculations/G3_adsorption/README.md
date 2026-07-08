# G3_adsorption layout

- `S1_Pd100/`, `S2_PdO101_Pd100/`, `S3_PdO100/`, `S3b_PdO100_PdOterm/`, `S4_PdO2_110/` —
  per-surface AutoAdsorbate candidate pools + MLIP (MACE) screening (`MLIP_phase1/2/3`,
  `*_filtered` = post geometry-refilter pass).
- `DFT_shortlist_v1/`, `DFT_shortlist_v2/`, `DFT_shortlist_v3/` — successive DFT-candidate
  picks from the MLIP-ranked pool (v1: original heuristic pick, v2: guide-strict repick,
  v3: energy + xy-distance dedup, **current/active** — this is what
  `scripts/setup_dft_v3.py` reads to build `calculations/T1_16_DFT_L2/`).
  Kept side by side for traceability; only v3 feeds live DFT jobs.
- `convB_binding/` — Convention-B (frozen-slab) binding-energy JSONs.
- `MLIP_phase*_summary.json`, `MLIP_phase*.log` — global (all-surface) MLIP run logs/summaries.
