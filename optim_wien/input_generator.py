"""
WIEN2k Input File Generator.

Generates valid case.in0, case.in1, case.in1c, case.in2, case.inm,
case.klist, case.in2c (hybrid functionals), and case.machines files
based on optimized parameters.

References:
  - WIEN2k User's Guide, Sections 4.2–4.6
  - WIEN2k FAQ: http://www.wien2k.at/reg_user/faq/
"""

import os


def generate_all_inputs(
    output_dir,
    basename,
    structure,
    rmt_result,
    rkmax_result,
    gmax_result,
    lmax_result,
    kmesh_result,
    mixing_result,
    core_valence_result,
    calc_type="scf",
    vxc_type=13,
    magnetic=False,
    spin_polarized=False,
    vxc_label=None,
    nproc_machines=0,
):
    os.makedirs(output_dir, exist_ok=True)
    files = {}

    in0 = _gen_in0(basename, gmax_result, vxc_type, spin_polarized, calc_type)
    fpath = os.path.join(output_dir, f"{basename}.in0")
    with open(fpath, "w") as f:
        f.write(in0)
    files["in0"] = fpath

    in1 = _gen_in1(basename, structure, rmt_result, rkmax_result,
                    lmax_result, core_valence_result, calc_type)
    fpath = os.path.join(output_dir, f"{basename}.in1")
    with open(fpath, "w") as f:
        f.write(in1)
    files["in1"] = fpath

    in1c = _gen_in1c(core_valence_result)
    fpath = os.path.join(output_dir, f"{basename}.in1c")
    with open(fpath, "w") as f:
        f.write(in1c)
    files["in1c"] = fpath

    in2 = _gen_in2(basename, mixing_result, magnetic, spin_polarized)
    fpath = os.path.join(output_dir, f"{basename}.in2")
    with open(fpath, "w") as f:
        f.write(in2)
    files["in2"] = fpath

    inm = _gen_inm(basename, mixing_result, magnetic, spin_polarized)
    fpath = os.path.join(output_dir, f"{basename}.inm")
    with open(fpath, "w") as f:
        f.write(inm)
    files["inm"] = fpath

    klist = _gen_klist(basename, kmesh_result)
    fpath = os.path.join(output_dir, f"{basename}.klist")
    with open(fpath, "w") as f:
        f.write(klist)
    files["klist"] = fpath

    if vxc_label and vxc_label.upper() == "HSE":
        in2c = _gen_in2c(basename, rkmax_result.rkmax, lmax_result)
        fpath = os.path.join(output_dir, f"{basename}.in2c")
        with open(fpath, "w") as f:
            f.write(in2c)
        files["in2c"] = fpath

    if nproc_machines > 0:
        from .submit import generate_machines
        mpath = generate_machines(basename, nproc_machines, output_dir)
        files["machines"] = mpath

    return files


def _gen_in0(case, gmax_result, vxc_type, spin_polarized, calc_type):
    mode = "FOR" if calc_type in ("forces",) else "TOT"
    lines = [
        f"{mode}             ({mode}/FOR with {case})",
        f"NR2V      IFFT      (R2V)",
        f"  {gmax_result.gmax:5.1f}       {gmax_result.gmax:5.1f}     "
        f"(GMAX for POTENTIAL and CHARGE)",
        "SPIN" if spin_polarized else "NON-SPIN",
        f"{vxc_type}  0.0    (VXCTYPE REL)",
        "",
    ]
    return "\n".join(lines)


def _gen_in1(case, structure, rmt_result, rkmax_result, lmax_result,
             core_valence_result, calc_type):
    ecut_abs = abs(core_valence_result.ecut)
    vnmt = _v_nmt(len(structure.atoms), core_valence_result.use_hdlo)
    global_lmax = max(lmax_result.lmax_values) if lmax_result.lmax_values else 6

    lines = [
        f"WFFIL        (WFPRI, SUPWF)",
        f" {rkmax_result.rkmax:5.1f}  {global_lmax:4d}  {vnmt:4d}  "
        f"(R-MT*K-MAX; MAX L IN WF; V-NMT)",
        f"  0.30    5  0      (GLOBAL E-PARAMETER with n other choices, "
        f"global APW/LAPW)",
        f" 0  -{ecut_abs:.1f}      0.0001      0.000     "
        f"EMIN  (semi-core electrons if negative)",
        f"           stop switch when EMIN reached",
        f"K-VECTORS FOR UNIT CELL WITH   {case}",
    ]

    for i, atom in enumerate(structure.atoms):
        rmt = rmt_result.rmt_values[i]
        lmax = lmax_result.lmax_values[i] if i < len(lmax_result.lmax_values) else 6
        lines.append(
            f" {lmax} 0 0    lmax for {atom.element}, atom {i+1}"
        )
        lines.append(
            f" 0.30    5  0     global parameters for {atom.element}"
        )
        lines.append(
            f" 1  -{ecut_abs:.1f}   0.00000   0.0000   "
            f"{atom.element} E-bottom {i+1}"
        )
        if core_valence_result.use_hdlo and rmt > 2.2:
            lines.append(
                f"LOCAL ORBITAL:    {atom.element}  ({i+1})   HDLO"
            )

    lines.append("")
    return "\n".join(lines)


def _v_nmt(num_atoms, use_hdlo):
    if use_hdlo or num_atoms >= 10:
        return 10
    if num_atoms >= 5:
        return 8
    return 6


def _gen_in1c(core_valence_result):
    ecut_abs = abs(core_valence_result.ecut)
    return "\n".join([
        "WFFIL        (WFPRI, SUPWF)",
        f" -{ecut_abs:.1f}  2    E-param for core states",
        "",
    ])


def _gen_in2(case, mixing_result, magnetic, spin_polarized):
    lines = [f"TOT             ({case})"]

    if mixing_result.scheme == "MSR1a":
        lines.append(
            f"MSR1a {mixing_result.mixing_factor:.2f}    "
            f"MIXING FACTOR FOR DENSITY"
        )
        lines.append(" 0.005  1.0   ")
        lines.append(" 0.0001  0.0001")
        lines.append("YES     (STORE restart file)")
    elif mixing_result.scheme == "MSEC1":
        lines.append(
            f"MSEC1 {mixing_result.mixing_factor:.2f}    "
            f"MIXING FACTOR FOR DENSITY"
        )
        lines.append(" 0.005  1.0   ")
        lines.append(" 0.0001  0.0001")
        lines.append("YES     (STORE restart file)")
    else:
        lines.append(
            f"PRATT  {mixing_result.mixing_factor:.2f}    "
            f"MIXING FACTOR FOR DENSITY"
        )

    if mixing_result.tetra:
        lines.append("TETRA    101")
    else:
        lines.append("GAUSS    0.0020    (GAUSS method)")

    lines.append(f" {mixing_result.temp:.4f}        TEMP")
    lines.append("")
    lines.append(f" {mixing_result.econv:.6f}    ENERGY CONVERGENCE")
    lines.append(f" {mixing_result.cconv:.6f}    CHARGE CONVERGENCE")

    if magnetic and spin_polarized:
        lines.append(f" {mixing_result.econv:.6f}    SPIN CONVERGENCE")
        lines.append(f" {mixing_result.cconv:.6f}    SPIN CHARGE CONVERGENCE")

    lines.append(f" {mixing_result.max_scf}      MAX SCF CYCLES")
    lines.append("")
    return "\n".join(lines)


def _gen_inm(case, mixing_result, magnetic, spin_polarized):
    if magnetic and spin_polarized:
        return "\n".join([
            f"{mixing_result.scheme} {mixing_result.mixing_factor:.2f}    "
            f"({case})",
            f"0.0 0.0 0.0    Bext",
            "",
        ])
    return f"NONE             ({case})\n"


def _gen_klist(case, kmesh_result):
    n1, n2, n3 = kmesh_result.mesh
    inv_flag = 1 if kmesh_result.add_inversion else 0
    s1, s2, s3 = kmesh_result.shift
    lines = [
        f" {case}",
        f" {n1:4d} {n2:4d} {n3:4d}    {kmesh_result.total_points}  "
        f"number of k-points",
        f"   -6  {inv_flag}     add INV" if inv_flag else f"   -6  0     no INV",
        f" {s1:5.2f} {s2:5.2f} {s3:5.2f}     shift",
    ]
    # Generate actual Monkhorst-Pack coordinates
    # k_i = (2*j - n_i + 1) / (2*n_i) + s_i / n_i  for j = 0..n_i-1
    weight = 1.0 / (n1 * n2 * n3)
    for j in range(n1):
        kx = (2.0 * j - n1 + 1) / (2.0 * n1) + s1 / n1
        for k in range(n2):
            ky = (2.0 * k - n2 + 1) / (2.0 * n2) + s2 / n2
            for l in range(n3):
                kz = (2.0 * l - n3 + 1) / (2.0 * n3) + s3 / n3
                lines.append(
                    f"{kx:12.8f} {ky:12.8f} {kz:12.8f} {weight:12.8f}"
                )
    lines.append("END")
    lines.append("")
    return "\n".join(lines)


def _gen_in2c(case, rkmax, lmax_result):
    """Generate case.in2c for hybrid functional (HSE) calculations.

    Hybrid functionals in WIEN2k require in2c for the exact-exchange
    part of the calculation. The FOCK operator replaces the Coulomb
    potential and requires a separate cutoff specification.

    Reference: WIEN2k User's Guide, Section 4.10 (hybrid functionals)
    """
    global_lmax = max(lmax_result.lmax_values) if lmax_result.lmax_values else 6
    lines = [
        f"TOT             ({case})",
        f" 0.30    5  0      (E-param with n other choices, global APW/LAPW)",
        f" {rkmax:5.1f}  {global_lmax:4d}        (RKMAX; GLOBAL LMAX FOR FOCK OPERATOR)",
        "  0    0    0        (FOCK-CONTROL: 0=normal Fock exchange)",
        "  0.25               (MIXING parameter for hybrid functional)",
        "",
    ]
    return "\n".join(lines)
