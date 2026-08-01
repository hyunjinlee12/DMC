# v5add16 bundle — verification report

Built 16 new candidate dirs under `calculations/T1_16_DFT_L2_v5add16/`

## Per-candidate verification

| # | sid | ads | idx | priority | POSCAR species | POTCAR TITEL | bash -n |
|---|---|---|---|---|---|---|---|
| 1 | S1 | coads | 487 | RECOMMENDED | `C H O Pd` | `C H O Pd` ✅ | OK |
| 2 | S1 | coads | 4021 | RECOMMENDED | `C H O Pd` | `C H O Pd` ✅ | OK |
| 3 | S2 | coads | 3633 | RECOMMENDED | `C H O Pd` | `C H O Pd` ✅ | OK |
| 4 | S2 | coads | 8079 | RECOMMENDED | `C H O Pd` | `C H O Pd` ✅ | OK |
| 5 | S3 | coads | 3481 | MUST | `C H O Pd` | `C H O Pd` ✅ | OK |
| 6 | S3 | coads | 5161 | MUST | `C H O Pd` | `C H O Pd` ✅ | OK |
| 7 | S3b | coads | 2051 | MUST-diagnostic | `C H O Pd` | `C H O Pd` ✅ | OK |
| 8 | S3b | coads | 2754 | MUST | `C H O Pd` | `C H O Pd` ✅ | OK |
| 9 | S1 | CO | 43 | OPTIONAL | `C O Pd` | `C O Pd` ✅ | OK |
| 10 | S1 | CH3O | 283 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |
| 11 | S2 | CH3O | 217 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |
| 12 | S2 | CH3O | 496 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |
| 13 | S3 | CH3O | 395 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |
| 14 | S3b | CO | 6 | OPTIONAL | `C O Pd` | `C O Pd` ✅ | OK |
| 15 | S4 | CH3O | 329 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |
| 16 | S4 | CH3O | 383 | OPTIONAL | `C H O Pd` | `C H O Pd` ✅ | OK |

## Issues found: 0
- None. All POSCAR/POTCAR species orders match; all submit scripts pass `bash -n`.

## Guardrails observed

- `T1_16_DFT_L2/` untouched (source of INCAR/POTCAR templates, read-only).
- No jobs submitted. Each dir contains submit_vasp_gpu.sh identical to v4 baseline;
  H200 environment paths must be adjusted before use, same as pending65 bundle.
- `manifest.csv` distinguishes 70 v4-existing (5 DONE, 65 PENDING) vs 16 v5-new (all PENDING).

## Special note

- **S3b/coads idx=2051** flagged `MUST-diagnostic` in metadata (ΔE=0.554 eV).
  Only reactive-close candidate available in the S3b coads MLIP pool. DFT will
  determine whether it remains bound at this configuration or relaxes to a
  lower-E arrangement.
