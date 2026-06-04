"""INCAR templates for Pd/PdO/PdO2 DMC project.

Usage:
    from incar_templates import get_incar
    incar_str = get_incar('bulk_metal')  # Pd
    incar_str = get_incar('bulk_oxide')  # PdO, PdO2
    incar_str = get_incar('slab_metal')
    incar_str = get_incar('slab_oxide')
"""

COMMON = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "LASPH": ".TRUE.",
    "ADDGRID": ".TRUE.",
    "ISPIN": 2,
    "IVDW": 12,
    "EDIFF": 1e-6,
    "NELM": 200,
    "NELMIN": 5,
    "ALGO": "Normal",
    "NCORE": 1,
    "LREAL": "Auto",
    "LWAVE": ".FALSE.",
    "LCHARG": ".FALSE.",
    "LORBIT": 11,
    "ISYM": 0,
}

TEMPLATES = {
    "bulk_metal": {
        **COMMON,
        "IBRION": 2,
        "NSW": 200,
        "ISIF": 3,
        "ISMEAR": 1,
        "SIGMA": 0.10,
        "EDIFFG": -0.01,
        "ISYM": 2,
        "LREAL": ".FALSE.",
    },
    "bulk_oxide": {
        **COMMON,
        "IBRION": 2,
        "NSW": 200,
        "ISIF": 3,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "EDIFFG": -0.01,
        "ISYM": 2,
        "LREAL": ".FALSE.",
    },
    "slab_metal": {
        **COMMON,
        "IBRION": 2,
        "NSW": 300,
        "ISIF": 2,
        "ISMEAR": 1,
        "SIGMA": 0.10,
        "EDIFFG": -0.03,
        "LDIPOL": ".TRUE.",
        "IDIPOL": 3,
    },
    "slab_oxide": {
        **COMMON,
        "IBRION": 2,
        "NSW": 300,
        "ISIF": 2,
        "ISMEAR": 0,
        "SIGMA": 0.05,
        "EDIFFG": -0.03,
        "LDIPOL": ".TRUE.",
        "IDIPOL": 3,
    },
    "solvation": {
        "LSOL": ".TRUE.",
        "EB_K": 32.6,
        "TAU": 0,
    },
    "neb": {
        "IMAGES": 5,
        "SPRING": -5,
        "LCLIMB": ".TRUE.",
        "IBRION": 3,
        "POTIM": 0,
        "NSW": 300,
        "EDIFFG": -0.05,
    },
}


def get_incar(template_name, overrides=None):
    params = dict(TEMPLATES[template_name])
    if overrides:
        params.update(overrides)
    lines = ["SYSTEM = pddmc"]
    for k, v in params.items():
        lines.append(f"{k} = {v}")
    return "\n".join(lines) + "\n"


def write_incar(path, template_name, overrides=None):
    with open(path, "w") as f:
        f.write(get_incar(template_name, overrides))
