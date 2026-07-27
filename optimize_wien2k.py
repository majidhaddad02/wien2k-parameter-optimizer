#!/usr/bin/env python3
"""
================================================================================
WIEN2k Comprehensive Parameter Optimization & Convergence Tool
================================================================================

Automatically optimizes ALL WIEN2k preprocessing parameters from a case.struct file:

  RMT     — Muffin-Tin radii (4 strict conditions)
  RKMAX   — Plane-wave cutoff (Blaha's reference table)
  GMAX    — Fourier expansion for density/potential
  LMAX/LVNS — Angular momentum cutoffs
  k-mesh  — Brillouin zone sampling (adaptive)
  Mixing  — SCF convergence scheme (PRATT/MSR1a/MSEC1)
  TEMP    — Fermi smearing for metals
  Core/Valence — Ecut and HDLO recommendations

Generates:
  - Comprehensive optimization report
  - case.in0, case.in1, case.in2, case.inm, case.klist input files
  - Optimized case.struct file

Usage:
    python optimize_wien2k.py case.struct [options]

Options:
    --calc-type TYPE     scf|relaxation|optimization|eos|forces|efg (default: scf)
    --precision PREC     screening|coarse|medium|high|very_high (default: medium)
    --refinement REF     coarse|medium|fine|very_fine (default: medium)
    --system-type TYPE   metal_small|semiconductor|insulator|... (auto-detect)
    --vxc TYPE           pbe|lda|wc|pbesol|scan (default: pbe)
    --magnetic           Enable spin-polarized calculation
    --output DIR         Output directory (default: ./optim_results)
    --no-input-files     Skip generating WIEN2k input files
    --quiet              Minimal console output

References:
    P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
    WIEN2k User's Guide & FAQ
"""

import argparse
import os
import sys

from optim_wien.constants import CalcType, Precision, VXCTYPE_PBE
from optim_wien.struct_parser import parse_struct
from optim_wien.rmt import optimize_rmt
from optim_wien.rkmax import optimize_rkmax
from optim_wien.gmax import optimize_gmax
from optim_wien.lmax import optimize_lmax
from optim_wien.kmesh import optimize_kmesh
from optim_wien.mixing import optimize_mixing
from optim_wien.core_valence import optimize_core_valence
from optim_wien.input_generator import generate_all_inputs
from optim_wien.report import generate_report, write_optimized_struct
from optim_wien.convergence import auto_converge, wien2k_available, ConvergenceResult


VXC_MAP = {
    "pbe": 13, "lda": 5, "wc": 11, "pbesol": 19, "scan": 28, "hse": 40,
}

CALC_MAP = {
    "scf": CalcType.SCF, "relaxation": CalcType.RELAXATION,
    "optimization": CalcType.OPTIMIZATION, "eos": CalcType.EOS,
    "forces": CalcType.FORCES, "efg": CalcType.EFG,
}

PREC_MAP = {
    "screening": Precision.SCREENING, "coarse": Precision.COARSE,
    "medium": Precision.MEDIUM, "high": Precision.HIGH,
    "very_high": Precision.VERY_HIGH,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="WIEN2k Comprehensive Parameter Optimization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python optimize_wien2k.py BaTiO3.struct
  python optimize_wien2k.py BaTiO3.struct --calc-type relaxation --precision high
  python optimize_wien2k.py Fe.struct --magnetic --precision very_high
  python optimize_wien2k.py slab.struct --system-type surface
        """,
    )
    parser.add_argument("struct_file", help="Path to case.struct file")
    parser.add_argument("--calc-type", default="scf",
                        choices=list(CALC_MAP.keys()), help="Calculation type")
    parser.add_argument("--precision", default="medium",
                        choices=list(PREC_MAP.keys()), help="Precision level")
    parser.add_argument("--refinement", default="medium",
                        choices=["coarse", "medium", "fine", "very_fine"],
                        help="k-mesh refinement")
    parser.add_argument("--system-type", default=None,
                        choices=["metal_small", "semiconductor", "insulator",
                                 "metal_large", "insulator_large", "surface",
                                 "molecule"])
    parser.add_argument("--vxc", default="pbe",
                        choices=list(VXC_MAP.keys()), help="XC functional")
    parser.add_argument("--magnetic", action="store_true",
                        help="Enable spin-polarized calculation")
    parser.add_argument("--auto-converge", action="store_true",
                        help="Auto-converge k-mesh & RKMAX by running WIEN2k")
    parser.add_argument("--output", default="./optim_results",
                        help="Output directory")
    parser.add_argument("--no-input-files", action="store_true",
                        help="Skip input file generation")
    parser.add_argument("--quiet", action="store_true",
                        help="Minimal output")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.struct_file):
        print(f"ERROR: {args.struct_file} not found.", file=sys.stderr)
        sys.exit(1)

    calc_type = CALC_MAP[args.calc_type]
    precision = PREC_MAP[args.precision]
    vxc = VXC_MAP.get(args.vxc, VXCTYPE_PBE)

    structure = parse_struct(args.struct_file)
    if not structure.atoms:
        print("ERROR: No atoms found.", file=sys.stderr)
        sys.exit(1)

    if not args.quiet:
        print(f"Structure: {structure.title}")
        print(f"  Lattice: {structure.lattice_type}")
        print(f"  a={structure.a:.4f} b={structure.b:.4f} c={structure.c:.4f}")
        print(f"  α={structure.alpha:.1f} β={structure.beta:.1f} γ={structure.gamma:.1f}")
        print(f"  Volume: {structure.volume:.2f} bohr³")
        print(f"  Atoms: {structure.num_atoms_primitive} "
              f"({len(structure.atoms)} non-equiv)")
        for a in structure.atoms:
            print(f"    {a.element}: Z={a.z}, mult={a.mult}, RMT={a.rmt:.3f}")
        print()

    if not args.quiet:
        print("Optimizing RMT...")
    rmt_result = optimize_rmt(structure, calc_type, precision)

    if not args.quiet:
        print("Optimizing RKMAX...")
    rkmax_result = optimize_rkmax(structure.atoms, rmt_result.rmt_values, precision)

    if not args.quiet:
        print("Optimizing GMAX...")
    gmax_result = optimize_gmax(structure.atoms, rmt_result.rmt_values, precision)

    if not args.quiet:
        print("Optimizing LMAX/LVNS...")
    lmax_result = optimize_lmax(structure.atoms, rmt_result.rmt_values)

    if not args.quiet:
        print("Optimizing k-mesh...")
    kmesh_result = optimize_kmesh(structure, refinement=args.refinement,
                                   system_type=args.system_type)

    if not args.quiet:
        print("Optimizing mixing/TEMP...")
    mixing_result = optimize_mixing(structure,
                                     system_type=kmesh_result.system_type,
                                     calc_type=calc_type, precision=precision,
                                     magnetic=args.magnetic)

    if not args.quiet:
        print("Optimizing core/valence...")
    core_valence_result = optimize_core_valence(structure.atoms,
                                                 rmt_result.rmt_values, precision)

    os.makedirs(args.output, exist_ok=True)
    basename = os.path.splitext(os.path.basename(args.struct_file))[0]
    if basename.endswith(".struct"):
        basename = basename[:-7]

    if not args.quiet:
        print("Generating report...")
    report = generate_report(
        structure, rmt_result, rkmax_result, gmax_result, lmax_result,
        kmesh_result, mixing_result, core_valence_result,
        calc_type=args.calc_type, precision=args.precision,
        struct_path=args.struct_file,
    )

    rpt_path = os.path.join(args.output, f"{basename}_optimization_report.txt")
    with open(rpt_path, "w") as f:
        f.write(report)

    struct_out = os.path.join(args.output, f"{basename}.struct_optimized")
    write_optimized_struct(structure, rmt_result.rmt_values, struct_out)

    gen_files = {}
    if not args.no_input_files:
        if not args.quiet:
            print("Generating WIEN2k input files...")
        gen_files = generate_all_inputs(
            args.output, basename, structure,
            rmt_result, rkmax_result, gmax_result, lmax_result,
            kmesh_result, mixing_result, core_valence_result,
            calc_type=args.calc_type, vxc_type=vxc,
            magnetic=args.magnetic, spin_polarized=args.magnetic,
        )

    conv_result = ConvergenceResult()
    if args.auto_converge:
        if not wien2k_available():
            if not args.quiet:
                print()
                print("⚠  WIEN2k not found in PATH. Cannot auto-converge.")
                print("   Install WIEN2k and source w2web environment first.")
                print("   Falling back to recommended parameters only.")
            conv_result = ConvergenceResult(
                warnings=["WIEN2k not available — using recommended parameters."]
            )
        else:
            if not args.quiet:
                print()
                print("=" * 60)
                print("AUTO-CONVERGENCE (Blaha hierarchy)")
                print("=" * 60)
                print("k-mesh convergence → RKMAX convergence")
                print(f"Threshold: ΔE < 0.1 mRy")
                print()

            conv_case_dir = os.path.abspath(args.output)
            conv_result = auto_converge(
                case_dir=conv_case_dir,
                basename=basename,
                structure=structure,
                initial_rmt_result=rmt_result,
                initial_rkmax_result=rkmax_result,
                initial_kmesh_result=kmesh_result,
                initial_mixing_result=mixing_result,
                initial_core_valence_result=core_valence_result,
                initial_gmax_result=gmax_result,
                initial_lmax_result=lmax_result,
                vxc_type=vxc,
                magnetic=args.magnetic,
                kmesh_threshold=0.0001,
                rkmax_threshold=0.0001,
            )

            if not args.quiet:
                print()
                print("=" * 60)
                print("CONVERGENCE RESULTS")
                print("=" * 60)
                print(f"  Converged: {conv_result.converged}")
                if conv_result.final_rmts:
                    print(f"  Final RMTs: {', '.join(f'{structure.atoms[i].element}={conv_result.final_rmts[i]:.4f}' for i in range(len(conv_result.final_rmts)))}")
                if conv_result.final_kmesh[0] > 0:
                    print(f"  Final k-mesh: {conv_result.final_kmesh[0]}×"
                          f"{conv_result.final_kmesh[1]}×"
                          f"{conv_result.final_kmesh[2]}")
                print(f"  Final RKMAX: {conv_result.final_rkmax}")
                print(f"  Runtime: {conv_result.total_runtime:.1f}s")
                if conv_result.rmt_history:
                    print("  RMT history:")
                    for entry in conv_result.rmt_history:
                        it = entry["iteration"]
                        rmts = entry["rmts"]
                        leak = entry.get("max_leak", "?")
                        rmt_str = ", ".join(
                            f"{structure.atoms[i].element}={rmts[i]:.4f}"
                            for i in range(len(rmts))
                        )
                        print(f"    iter {it}: RMT=({rmt_str}), max_leak={leak}")
                        if "action" in entry:
                            print(f"      → {entry['action']}")
                if conv_result.kmesh_history:
                    print("  k-mesh history:")
                    for entry in conv_result.kmesh_history:
                        m = entry["mesh"]
                        e = entry["energy"]
                        es = f"{e:.8f}" if e else "N/A"
                        print(f"    {m[0]}×{m[1]}×{m[2]} → E = {es}")
                if conv_result.rkmax_history:
                    print("  RKMAX history:")
                    for entry in conv_result.rkmax_history:
                        r = entry["rkmax"]
                        e = entry["energy"]
                        es = f"{e:.8f}" if e else "N/A"
                        print(f"    RKMAX={r} → E = {es}")
                for w in conv_result.warnings:
                    print(f"  ⚠ {w}")
                print()

    if not args.quiet:
        print()
        print(report)
        print()
        print(f"Output directory: {args.output}/")
        print(f"  {basename}_optimization_report.txt")
        print(f"  {basename}.struct_optimized")
        if not args.no_input_files:
            for name, fpath in gen_files.items():
                print(f"  {basename}.{name} → {name}")
        print()

    has_critical = any("CRITICAL" in w for w in rmt_result.core_leakage_warnings)
    has_final = any("FINAL WARNING" in w for w in rmt_result.overlap_warnings)

    if has_critical or has_final:
        print("⚠  Critical warnings — manual review required.")
        sys.exit(2)

    if not args.quiet:
        n1, n2, n3 = kmesh_result.mesh
        print(f"Ready for WIEN2k initialization:")
        print(f"  init_lapw -b -rkmax {rkmax_result.rkmax} "
              f"-numk {n1*100+n2*10+n3} "
              f"-ecut {int(abs(core_valence_result.ecut))}")
        print(f"  run_lapw -p")
        print()


if __name__ == "__main__":
    main()
