# Judge Committee — quality control for the DMC project

A blind 3-judge committee + chair pattern (inspired by LLM-as-Judge and Mixture-of-Agents
literature) sits between **producers** and **orchestrator** to QC producer outputs
before they advance project gates.

## When to convene

Convene the committee whenever a producer agent finishes a high-stakes output
that will gate downstream work:

| Trigger | Producer | Convene? |
|---|---|---|
| MACE ranking phase done (Phase 1/2/3) | ml-trainer | **YES** — before picking DFT shortlist |
| DFT shortlist results in | simulation | **YES** — before T1.18 ads E table |
| NEB profile done | simulation | **YES** — before benchmark vs Shi 2024 |
| Slab built for a new termination | data-curator | optional (cheap retry) |
| Routine docs / refactor | github-manager | NO |
| Lit search | literature | NO |

## Procedure (orchestrator runs this)

1. Producer returns artifacts.
2. Dispatch the THREE judges **in parallel** (single message, three Agent calls),
   each given the same producer artifacts list — and ONLY artifacts, never each
   other's verdicts.
3. Collect three verdict JSONs.
4. Dispatch `committee-chair` with the three verdicts as input.
5. Read chair output. Apply decision:
   - **GO** → proceed; record asks to track during the next phase.
   - **REVIEW** → execute consolidated asks before proceeding (often: producer
     re-run on a specific subset; an additional sanity check; advisor consult).
   - **REJECT** → producer must fix critical violation. Re-dispatch with the
     specific corrective brief from chair's `critical_violations`.
6. Persist chair output under `calculations/<phase>/committee_<phase>_<date>.json`.

## Why blind judges
Judges seeing each other's verdicts collapse to groupthink. Independence is the
whole point — three different "biases" (methods, physics, statistics) approximate
an unbiased aggregate. If you only had one judge, you'd inherit its blind spots.

Once judges have spoken, the chair sees everything and can identify **disagreement**
— often the most informative signal (e.g., methods passes but physics rejects =
producer respected the recipe but the chemistry is wrong; investigate the
producer's INPUT, not the producer).

## Judge roles (summary — Paimon-extended)

- **`judge-methods`** — Did the producer follow project rules? (config, env,
  gate prerequisites, file hygiene). Reads CLAUDE.md, memory, workplan, then
  reads the run. Output: rule compliance + specific violations.
- **`judge-physics`** — Is the chemistry plausible? (bond lengths, broken
  adsorbates, exploded slabs, E in expected literature range). Opens the actual
  structures and spot-checks. Output: physical sanity + concerns.
- **`judge-statistics`** — Does the data distribution support the conclusion?
  (convergence rate, dedup ratio, energy spread, coverage of bins). Computes
  metrics from JSON / traj. Output: quantitative validity + suspect metrics.
- **`judge-silent-error`** (Paimon 2026 §2.2 inspired) — Catches "looks right
  but is wrong": claimed metrics that don't match recomputed values, sort by
  wrong key, broken units, mislabeled "DFT" vs "MLIP" energies. Recomputes
  producer's claims from raw artifacts.
- **`judge-malicious`** (Paimon §3.3 inspired) — Binary trust audit. Detects
  premature completion, fabricated numbers, bypassed checks. If fires
  → committee Fail (no retry, escalate).
- **`committee-chair`** — Aggregates 5 verdicts → 4-label decision
  (Pass / Pass-with-caveats / Revise / Fail) using priority logic. Hands
  Revise to summarizer. Surfaces disagreement patterns.

## Supporting agents (not judges)

- **`annotator`** — Pre-judge / pre-producer context enricher. Reads docs,
  highlights governing sections verbatim, makes implicit assumptions explicit.
  Use when starting a new T-number or before high-stakes judge convocation.
- **`summarizer`** — When chair returns Revise, condenses raw judge JSON into
  1-3 actionable items the producer can execute. Reconciles contradictions,
  prioritizes by severity.

## 4-label verdict (Paimon-style)

| Label | Meaning | Action |
|---|---|---|
| **Pass** | All judges Pass | Proceed |
| **Pass-with-caveats** | 1-2 Concerns, rest Pass | Proceed + log concerns to STATUS.md |
| **Revise** | Any Reject OR >50% Concerns | Summarizer → producer fix → re-judge (retry budget 1-2) |
| **Fail** | judge-malicious=Malicious | STOP, escalate to user/advisor |

## Revise loop (Paimon §3.3)

```
chair returns Revise
  ↓
summarizer condenses asks → 1-3 actionable items (prioritized)
  ↓
producer re-dispatched with the items as brief
  ↓
producer fixes + re-runs (or re-runs only the affected part if resumable)
  ↓
committee re-convenes (only the judges whose verdict was non-Pass)
  ↓
If still Revise after retry_budget exhausted → Fail
If Pass / Pass-with-caveats → proceed
```

## Producer-type weighting (chair's job)

Different producer outputs need different judge emphases:

```
ml-trainer (MLIP ranking)     methods 30%, physics 30%, statistics 40%
data-curator (structures)     methods 30%, physics 50%, statistics 20%
simulation (DFT/NEB)          methods 40%, physics 50%, statistics 10%
analyst (interpretation)      methods 20%, physics 50%, statistics 30%
```

Statistics gets less weight on a single-shot DFT relax (n=1, no distribution),
much more on a 38k MLIP ranking (the whole point is the distribution).

## Memory: when judges disagree

If `judge-methods` says GO and `judge-physics` says REJECT, that means:
- Producer followed the recipe ✓
- But the OUTPUT chemistry is broken ✗
→ Likely the producer's **INPUT** was bad (e.g., wrong structure handed by
data-curator, wrong reference numbers, lit ref outdated).
→ Orchestrator should re-examine upstream, not blame the producer.

If `judge-physics` says GO and `judge-methods` says REJECT:
- Output looks chemically real, but recipe was wrong.
→ Maybe right answer by luck — DON'T accept. Re-run with correct recipe.

## Cost discipline

Convening the committee = 4 extra agent invocations. Reserve for high-stakes
gate transitions, not every routine task. Typical pattern:

- Phase 1 done → committee (4 agents)
- Phase 1 results sound → Phase 2 dispatch (no committee yet)
- Phase 2 done → committee
- DFT shortlist DFT run → committee per surface result
- ...

If you committee everything you 5×-multiply your token budget.

## Example invocation (orchestrator pseudo-code)

```python
# After ml-trainer Phase 1 done:
artifacts = "calculations/G3_adsorption/MLIP_phase1_summary.json + per-surface MLIP_phase1/ dirs"
v1 = Agent("judge-methods",    prompt=f"Judge ml-trainer Phase 1 output: {artifacts}")
v2 = Agent("judge-physics",    prompt=f"Judge ml-trainer Phase 1 output: {artifacts}")  # in parallel
v3 = Agent("judge-statistics", prompt=f"Judge ml-trainer Phase 1 output: {artifacts}")  # in parallel

# Wait for all three
verdict = Agent("committee-chair", prompt=f"Aggregate: v1={v1}, v2={v2}, v3={v3}, producer_type=ml-trainer")

# Decide
if verdict.decision == "GO":
    dispatch Phase 2
elif verdict.decision == "REVIEW":
    address verdict.consolidated_asks first
else:
    rerun ml-trainer with verdict.critical_violations fixed
```
