"""After all 11 seed dirs finish, run this to summarize:
- final energy(sigma->0)
- final total magnetization
- per-orbital tot magnetization (s, p, d)
- electronic convergence status (for NSW=0 static)
- fatal error detection

For static (NSW=0) calcs, "reached required accuracy" (ionic loop marker) is
NOT expected. Electronic convergence in VASP prints
"aborting loop because EDIFF is reached" (or the equivalent for the ALGO used).

Usage:
    python analyze_seed_results.py > summary.txt
"""
import re
from pathlib import Path

for d in sorted(Path('magnetic_seed_test').iterdir()):
    if not d.is_dir(): continue
    outcar = d/'OUTCAR'
    if not outcar.exists():
        print(f'{d.name:<38} STATUS = not yet run')
        continue
    text = outcar.read_text(errors='replace')
    # electronic convergence marker (works for static and relax)
    ediff_reached = 'aborting loop because EDIFF is reached' in text
    # generic completion (also present at end of static)
    finished = ('General timing and accounting informations' in text
                or 'Total CPU time used' in text)
    # any fatal error
    fatal_patterns = [r'ERROR', r'BRMIX', r'Sub-Space-Matrix is not hermitian',
                      r'ZBRENT', r'WARNING: DENTET']
    fatal = any(re.search(p, text) for p in fatal_patterns)
    # extract energies
    E = re.findall(r'energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)', text)
    mag = re.findall(r'number of electron.*?magnetization\s+(-?\d+\.\d+)', text)
    E_final = float(E[-1]) if E else None
    m_final = float(mag[-1]) if mag else None
    status = ('OK' if (ediff_reached and finished and not fatal)
              else 'CHECK' if fatal
              else 'INCOMPLETE')
    print(f'{d.name:<38} status={status:<10} E={E_final} mu_tot={m_final} '
          f'EDIFF_reached={ediff_reached} finished={finished} fatal={fatal}')
