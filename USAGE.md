# USAGE — research-pd-dmc

Setup, how to launch this project, the gated workflow, how to talk to the Director,
running on a server, and cost control.

---

## 1. One-time setup

### 1a. Claude Code + GitHub
```bash
npm i -g @anthropic-ai/claude-code
claude --version
gh auth login                 # github-manager uses the gh CLI
```

### 1b. Materials + simulation Python environment
```bash
conda create -n pddmc python=3.11 -y && conda activate pddmc
pip install ase pymatgen mp-api          # structures, slab generation, MP retrieval
pip install rdkit autoadsorbate          # adsorbate generation
pip install mace-torch                   # MACE foundation MLIP (ranking) — GPU recommended
pip install numpy scipy pandas matplotlib
# VASP (+VASPsol) and VTST tools are provided by your cluster, not pip.
```

### 1c. Secrets via environment (never in files)
```bash
export MP_API_KEY="...your Materials Project key..."   # data-curator reads this
# add to ~/.bashrc / your job profile so sessions inherit it
```

### 1d. Cluster
The box where Claude Code runs must reach your scheduler (`sbatch`/`squeue`/`sacct`),
and VASP + VASPsol + POTCAR pseudopotentials + VTST scripts must be configured. The
`simulation` agent submits async jobs and polls them.

### 1e. Reference docs
Copy the four guideline files into `docs/` (see `docs/README.md`):
`DMC_Pd_workplan.md`, `DMC_Pd_package_guideline.md`,
`DMC_Pd_references_summary.md`, `DMC_Pd_student_guide.md`.
The Director and the agents read these as authoritative context.

---

## 2. Launch
```bash
cd research-pd-dmc
export CLAUDE_CODE_SUBAGENT_MODEL="claude-sonnet-4-5"   # workers on Sonnet
claude --model claude-opus-4-6                          # Director on Opus
```
The Director reads `PROJECT.md` + `docs/`, plans, dispatches specialists, and reports back.
(Newer model strings can be substituted if available.)

---

## 3. The gated workflow (don't skip gates)

This study advances through hard gates; the Director will pause at each for your review.

**Phase 1**
1. `literature` → scope the anchor refs + benchmarks + surface-model justification.
2. `data-curator` → fetch fcc Pd / tetragonal PdO / rutile-like PdO₂ bulk.
3. `simulation` → bulk relax + ENCUT/k-mesh convergence; lattice vs experiment. **→ G1 (checkpoint A)**
4. `data-curator` → build S1–S4 clean slabs from the relaxed bulk (asymmetric, vac 20 Å).
5. `simulation` → slab relax; `analyst` → rumpling/Pd–O/dipole/termination validation. **→ G2 (checkpoint B)**
6. `data-curator` → site maps + CO*/CH₃O* (`mode='all'`) + co-adsorption (site-pair wrapper).
7. `ml-trainer` → MACE relax + rank → per-surface DFT shortlist (distance-bin reps).
8. `simulation` → DFT adsorption: Level 1 vacuum → Level 2 VASPsol.
9. `analyst` → adsorption table + descriptor map + Case A–D → pick Phase-2 surfaces. **→ G3 (checkpoint E)**

**Phase 2** (selected surfaces only)
10. `data-curator` → reactive pairs + intermediates (CH₃OCO*, DMC*).
11. `ml-trainer` → rank → endpoint shortlist; `simulation` → endpoint opt → TS1/TS2 (CI-NEB) + side-path.
12. `analyst` → DMC Gibbs profile + Shi-2024 benchmark + side-path competition → conclusion. **→ G4**
13. `github-manager` → commit artifacts + open PR with the numbers.

---

## 4. Talking to the Director — patterns that work
- **One gated stage:** *"We're at G1. Have data-curator fetch the three bulk structures,
  then simulation relax them + run convergence; report lattice errors. Stop at G1."*
- **Call a specialist directly:** *"Use data-curator to generate CO* candidates on the
  S2 slab with mode='all', overlap_thr=1.25, and report the site distribution."*
- **Screen then judge:** *"Have ml-trainer MACE-rank the S1 co-adsorption pool and propose
  a 5-structure DFT shortlist spanning distance bins; then simulation runs Level 1 on them."*
- **Async DFT:** *"Submit the S3 slab relaxations; give me job IDs, we'll check back."*
  Later: *"Poll those jobs; when converged, send energies to analyst for the descriptor map."*
- **Stop points (cost):** *"Don't submit any NEB without showing me the cores×walltime estimate first."*

Remember: MLIP output is screening only — insist that conclusions rest on DFT.

---

## 5. Running on a server (remote control + persistence)
```bash
tmux new -s director
cd research-pd-dmc && export CLAUDE_CODE_SUBAGENT_MODEL="claude-sonnet-4-5"
claude --model claude-opus-4-6
# inside the session:
/rc          # Remote Control: drive from phone / claude.ai/code
/mobile      # QR to install the app
# detach: Ctrl-b d   reattach: tmux attach -t director
```
Code stays on the server; `/rc` is an outbound control channel. For sturdiness, wrap the
launch in a systemd user service so it restarts on reboot.

---

## 6. Cost control
- **Model split:** Director on Opus (planning/judgment), workers on Sonnet via
  `CLAUDE_CODE_SUBAGENT_MODEL`. `literature` can run on a cheaper model.
- **DFT is the real cost, not tokens.** Always get `simulation`'s estimate before batches;
  keep the MLIP funnel tight (2–4 DFT structures per surface/species) so DFT stays small.
- **Async simulation** keeps the Director from idling (billing) while jobs run on the cluster.
- **Parallel = linear cost:** only parallelize genuinely independent work (bulk T1.1–1.3;
  per-surface slabs; CO* ∥ CH₃O*; TS1 ∥ TS2 ∥ side-path).
