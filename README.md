# research-pd-dmc — multi-agent DFT study (Pd/PdO/PdO₂ DMC formation)

A terminal-based multi-agent Claude Code setup, adapted from `materials-ml-kit` for a
**DFT mechanistic catalysis study** (not ML property prediction). Question: does anodic
oxidation of Pd {100} (Pd⁰ → Pd²⁺ → Pd⁴⁺) promote or suppress electrochemical DMC
formation? See **`PROJECT.md`** for the brief and **`USAGE.md`** for the walkthrough.

## Current status
- **G1 (bulk)** ✅ passed · **G2 (slab)** ✅ passed — 5 surfaces validated (S1 Pd(100),
  S2 PdO(101)/Pd(100), S3 PdO(100), S3b PdO(100) PdO-term, S4 PdO2(110)).
- **G3 (adsorption), in progress**: MLIP (MACE) screening done on all 5 surfaces ×
  3 adsorbate types (CO*, CH₃O*, co-ads) → top-5 low-energy & site-diverse candidates
  per surface/adsorbate shortlisted (70 total, `T1_16_DFT_L2`).
  DFT (VASPsol) validation of the shortlist is underway — S1 Pd(100) CO* done (5/5);
  remaining 65 candidates queued.
- G3 gate (descriptor map + Case A–D) not yet reached.

## The team (strict hierarchy — everything routes through the Director)

```
              you ──orders──▶ DIRECTOR ──reports──▶ you
                                 │
 ┌──────────┬──────────────┬─────┴─────┬───────────┬────────────────┐
literature  data-curator   simulation  ml-trainer   analyst    github-manager
(refs +     (structures +  (VASP DFT:  (MACE MLIP   (descriptor (save / PR)
 benchmarks) slabs +        bulk/slab/  ranking —    map, Gibbs
             AutoAdsorbate) adsorption/ screening    profile,
                            NEB·async)  only)        Case A–D)
```

Pipeline is **gated** (G1–G4), not a single linear pass:

```
literature → data-curator(bulk) → simulation(bulk relax) ─[G1]─
→ data-curator(slabs) → simulation(slab relax)+analyst(validate) ─[G2]─
→ data-curator(adsorbates, mode='all') → ml-trainer(MACE rank) → simulation(DFT L1→L2 VASPsol)
→ analyst(descriptor map, Case A–D) ─[G3]─
→ [Phase 2 on selected surfaces] data-curator → ml-trainer → simulation(endpoints, TS1/TS2, side-path)
→ analyst(Gibbs profile, Shi-2024 benchmark) ─[G4]─ → github-manager(save + PR)
```

## Layout
```
research-pd-dmc/
├── CLAUDE.md                 # Director (reads PROJECT.md + docs/)
├── PROJECT.md                # the brief — edit this to retune the project
├── README.md
├── USAGE.md                  # detailed walkthrough
├── docs/                     # drop the 4 source guideline .md files here (see docs/README.md)
└── .claude/agents/
    ├── literature.md
    ├── data-curator.md       # structure & adsorbate builder
    ├── simulation.md         # VASP(+VASPsol) DFT engine
    ├── ml-trainer.md         # MACE MLIP screener (ranking only)
    ├── analyst.md            # descriptor map / Gibbs profile interpreter
    └── github-manager.md
```

## Quickstart
```bash
# one-time
npm i -g @anthropic-ai/claude-code
gh auth login
export MP_API_KEY="...your Materials Project key..."

# copy your four guideline docs into docs/ (see docs/README.md), then:
cd research-pd-dmc
export CLAUDE_CODE_SUBAGENT_MODEL="claude-sonnet-4-5"
claude --model claude-opus-4-6
```
Then talk to the Director in plain language, e.g. *"We're at T1.1. Fetch the three bulk
structures and run bulk relaxation + convergence; report lattice vs experiment for G1."*

> Note: the agent files here are tuned for DFT catalysis. To reuse them for other
> catalysis projects, copy `.claude/agents/*` back into your `materials-ml-kit/template/`.
