#!/usr/bin/env python3
"""
WIEN2k Comprehensive Parameter Optimization & Convergence Tool.

Two modes:
  1. Interactive wizard:  python optimize_wien2k.py -i
  2. Command-line:         python optimize_wien2k.py case.struct [options]

References:
    P. Blaha et al., J. Chem. Phys. 152, 074101 (2020)
    WIEN2k User's Guide & FAQ
"""

import argparse
import os
import sys
import time

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
from optim_wien.converge import (
    build_convergence_report, ConvergenceEngineResult, DEFAULT_ETOL,
)
from optim_wien.cli import (
    InteractiveWizard, ProgressTracker,
    box, header, section, info, warn, error, success, style, clear, echo,
    confirm, menu, banner, term_width, BOX,
)

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
  python optimize_wien2k.py -i                          # interactive wizard
  python optimize_wien2k.py BaTiO3.struct                # quick run
  python optimize_wien2k.py BaTiO3.struct --precision high --refinement fine
  python optimize_wien2k.py Fe.struct --magnetic --precision very_high
  python optimize_wien2k.py case.struct --converge rkmax,kmesh  # verify convergence
  python optimize_wien2k.py case.struct --auto-converge  # converge with WIEN2k
        """,
    )
    parser.add_argument(
        "struct_file", nargs="?", help="Path to case.struct file"
    )
    parser.add_argument(
        "-i", "--interactive", action="store_true",
        help="Launch interactive wizard mode"
    )
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
    parser.add_argument("--bandgap", type=float, default=None,
                        metavar="eV",
                        help="Known band gap (eV) — affects k-mesh density "
                             "(semiconductor=1.0, insulator=0.5, metal=0.0)")
    parser.add_argument("--vxc", default="pbe",
                        choices=list(VXC_MAP.keys()), help="XC functional")
    parser.add_argument("--magnetic", action="store_true",
                        help="Enable spin-polarized calculation")
    parser.add_argument("--auto-converge", action="store_true",
                        help="Auto-converge k-mesh & RKMAX by running WIEN2k")
    parser.add_argument("--converge", default=None, nargs="*",
                        choices=["rmt", "rkmax", "kmesh", "gmax"],
                        metavar="PARAM",
                        help="Convergence-verified optimization. "
                             "Comma-separated list of parameters: "
                             "rmt rkmax kmesh gmax (default: rkmax kmesh)")
    parser.add_argument("--etol", type=float, default=DEFAULT_ETOL,
                        metavar="mRy/atom",
                        help=f"Convergence tolerance in mRy/atom "
                             f"(default: {DEFAULT_ETOL})")
    parser.add_argument("--cluster-submit", action="store_true",
                        help="Submit convergence jobs to HPC scheduler "
                             "(requires submit callback)")
    parser.add_argument("--converge-report", default=None, metavar="PATH",
                        help="Path for convergence report (markdown). "
                             "Default: <output>/<case>_convergence_report.md")
    parser.add_argument("--output", default="./optim_results",
                        help="Output directory")
    parser.add_argument("--no-input-files", action="store_true",
                        help="Skip input file generation")
    parser.add_argument("--quiet", action="store_true",
                        help="Minimal output")
    parser.add_argument("--only", default=None, nargs="*",
                        choices=["rmt", "rkmax", "gmax", "lmax", "kmesh",
                                 "mixing", "core"],
                        metavar="STEP",
                        help="Run only specified step(s): rmt rkmax gmax "
                             "lmax kmesh mixing core")
    parser.add_argument("--strict-faq", action="store_true",
                        help="Strict FAQ mode: use Blaha table values directly, "
                             "no precision offsets. References: WIEN2k FAQ RMT, "
                             "RKMAX, kgen (http://www.wien2k.at/reg_user/faq/)")
    parser.add_argument("--prec", type=int, default=None,
                        choices=[0, 1, 2, 3],
                        help="init_lapw -prec flag value (0=fast, 1=standard, "
                             "2=high, 3=very high). Overrides --precision mapping. "
                             "Reference: WIEN2k User's Guide, Section 3.2")
    parser.add_argument("--machines", type=int, default=0, metavar="NPROC",
                        help="Generate case.machines file for NPROC parallel "
                             "processes. Reference: WIEN2k User's Guide, Section 3.5")
    parser.add_argument("--submit", default=None,
                        choices=["slurm", "pbs"], metavar="SCHEDULER",
                        help="Generate job submission script (SLURM or PBS) "
                             "with optimized parameters")
    return parser.parse_args()


def run_interactive():
    wizard = InteractiveWizard()
    config = wizard.run()
    if config is None:
        echo()
        echo(style("  Cancelled.", fg="bright_black"))
        return

    clear()
    banner()
    _run_optimization(config, interactive=True)


_STEP_DEPS = {
    "rmt":     frozenset(),
    "rkmax":   frozenset({"rmt"}),
    "gmax":    frozenset({"rmt"}),
    "lmax":    frozenset({"rmt"}),
    "kmesh":   frozenset(),
    "mixing":  frozenset({"kmesh"}),
    "core":    frozenset({"rmt"}),
}

_ALL_STEPS = ("rmt", "rkmax", "gmax", "lmax", "kmesh", "mixing", "core")


def _resolve_steps(only_flag):
    if only_flag is None or len(only_flag) == 0:
        return set(_ALL_STEPS)
    if isinstance(only_flag, str):
        requested = {only_flag}
    else:
        requested = set(only_flag)
    resolved = set(requested)
    for s in requested:
        deps = _STEP_DEPS.get(s, frozenset())
        resolved |= deps
    return resolved


def _run_optimization(config, interactive=False):
    args = type("Args", (), config)()

    active = _resolve_steps(getattr(args, "only", None))
    is_partial = active != set(_ALL_STEPS)

    if not os.path.isfile(args.struct_file):
        error(f"'{args.struct_file}' not found.")

    calc_type = CALC_MAP[args.calc_type]
    precision = PREC_MAP[args.precision]
    vxc = VXC_MAP.get(args.vxc, VXCTYPE_PBE)
    quiet = args.quiet
    verbose = not quiet

    structure = parse_struct(args.struct_file)
    if not structure.atoms:
        error("No atoms found in struct file.")

    strict_faq = getattr(args, "strict_faq", False)

    rmt_result = rkmax_result = gmax_result = lmax_result = None
    kmesh_result = mixing_result = core_valence_result = None

    _get_rmt = lambda: rmt_result.rmt_values  # noqa: E731
    _get_kmesh_stype = lambda: kmesh_result.system_type  # noqa: E731

    tracker = ProgressTracker()

    n_vis = len([s for s in ("rmt", "rkmax", "gmax", "lmax", "kmesh", "mixing", "core")
                 if s in active]) + 2  # +2 for Report and Input Files

    if verbose:
        clear()
        if is_partial:
            s_list = ", ".join(sorted(active))
            header(f"SELECTED STEPS: {s_list}", "bright_magenta")
        else:
            banner()

        header("STRUCTURE ANALYSIS", "bright_green")
        box(
            f"  {style('Title:', fg='bright_black', bold=True)}  "
            f"{style(structure.title, fg='white')}\n"
            f"  {style('Lattice:', fg='bright_black', bold=True)}  "
            f"{style(structure.lattice_type, fg='bright_magenta')}\n"
            f"  {style('a, b, c:', fg='bright_black', bold=True)}  "
            f"{style(f'{structure.a:.4f}  {structure.b:.4f}  {structure.c:.4f}  bohr', fg='white')}\n"
            f"  {style('α, β, γ:', fg='bright_black', bold=True)}  "
            f"{style(f'{structure.alpha:.1f}°  {structure.beta:.1f}°  {structure.gamma:.1f}°', fg='white')}\n"
            f"  {style('Volume:', fg='bright_black', bold=True)}  "
            f"{style(f'{structure.volume:.2f} bohr³', fg='bright_magenta')}\n"
            f"  {style('Atoms:', fg='bright_black', bold=True)}  "
            f"{style(f'{structure.num_atoms_primitive} ({len(structure.atoms)} non-equiv)', fg='white')}",
            fg="bright_green",
        )
        for a in structure.atoms:
            info(f"{a.element}", f"Z={a.z}  mult={a.mult}  initial RMT={a.rmt:.3f}")

        tracker.start(n_vis)

    t_start = time.time()

    def _step(name, fn, *, show=None):
        """Run one optimization step. If `show` is False, run silently (dependency)."""
        visible = show if show is not None else (name in active)
        if visible and verbose:
            tracker.step(name)
        result = fn()
        if visible and verbose:
            tracker.done()
        return result

    rmt_result = _step("rmt",
        lambda: optimize_rmt(structure, calc_type, precision))
    rkmax_result = _step("rkmax",
        lambda: optimize_rkmax(structure.atoms, _get_rmt(), precision, strict_faq=strict_faq))
    gmax_result = _step("gmax",
        lambda: optimize_gmax(structure.atoms, _get_rmt(), precision, strict_faq=strict_faq))
    lmax_result = _step("lmax",
        lambda: optimize_lmax(structure.atoms, _get_rmt()))
    kmesh_result = _step("kmesh",
        lambda: optimize_kmesh(structure,
                                refinement="medium" if strict_faq else args.refinement,
                                system_type=args.system_type,
                                bandgap=getattr(args, "bandgap", None)))
    mixing_result = _step("mixing",
        lambda: optimize_mixing(structure, system_type=_get_kmesh_stype(),
                                calc_type=calc_type, precision=precision,
                                magnetic=args.magnetic))
    core_valence_result = _step("core",
        lambda: optimize_core_valence(structure.atoms, _get_rmt(), precision))

    os.makedirs(args.output, exist_ok=True)
    basename = os.path.splitext(os.path.basename(args.struct_file))[0]
    if basename.endswith(".struct"):
        basename = basename[:-7]

    if verbose:
        tracker.step("report")
    report = generate_report(
        structure, rmt_result, rkmax_result, gmax_result, lmax_result,
        kmesh_result, mixing_result, core_valence_result,
        calc_type=args.calc_type, precision=args.precision,
        struct_path=args.struct_file,
    )
    if verbose:
        tracker.done()

    rpt_path = os.path.join(args.output, f"{basename}_optimization_report.txt")
    with open(rpt_path, "w") as f:
        f.write(report)

    struct_out = os.path.join(args.output, f"{basename}.struct_optimized")
    write_optimized_struct(structure, rmt_result.rmt_values, struct_out)

    gen_files = {}
    if not args.no_input_files:
        if verbose:
            tracker.step("input_files")
        gen_files = generate_all_inputs(
            args.output, basename, structure,
            rmt_result, rkmax_result, gmax_result, lmax_result,
            kmesh_result, mixing_result, core_valence_result,
            calc_type=args.calc_type, vxc_type=vxc,
            magnetic=args.magnetic, spin_polarized=args.magnetic,
            vxc_label=args.vxc,
            nproc_machines=getattr(args, "machines", 0),
        )
        if verbose:
            tracker.done()

    # ── Generate job submission scripts if requested ──
    gen_submit = getattr(args, "submit", None)
    if gen_submit:
        from optim_wien.submit import generate_slurm_script, generate_pbs_script
        n1, n2, n3 = kmesh_result.mesh
        nproc = max(getattr(args, "machines", 0), 4)
        kWargs = dict(
            basename=basename, work_dir=os.path.abspath(args.output),
            rkmax=rkmax_result.rkmax, numk=n1 * n2 * n3,
            ecut=int(abs(core_valence_result.ecut)), parallel=True,
            output_dir=args.output,
        )
        if gen_submit == "slurm":
            spath = generate_slurm_script(nproc=nproc, **kWargs)
            gen_files["slurm"] = spath
        elif gen_submit == "pbs":
            ppath = generate_pbs_script(nproc=nproc, **kWargs)
            gen_files["pbs"] = ppath

    if verbose:
        tracker.finish()

    conv_result = ConvergenceResult()
    conv_engine = None

    # ── NEW: Convergence-verified optimization (Aitken extrapolation) ──
    if args.converge is not None:
        converge_params_set = set(args.converge) if args.converge else {"rkmax", "kmesh"}
        for_forces = args.calc_type in ("forces", "relaxation")

        if not wien2k_available():
            if not quiet:
                warn("WIEN2k not found. Cannot run convergence engine.")
                echo(style("       Use --auto-converge with WIEN2k installed.", fg="bright_black"))
        else:
            if not quiet:
                header("CONVERGENCE-VERIFIED OPTIMIZATION", "bright_magenta")
                info("Method", "Aitken Δ² extrapolation with SCF confirmation")
                info("Parameters", ", ".join(sorted(converge_params_set)))
                info("Tolerance", f"{args.etol:.2f} mRy/atom")
                if for_forces:
                    etol_eff = args.etol / 10.0
                    info("Forces mode", f"etol tightened to {etol_eff:.2f} mRy/atom "
                         "(Blaha 2020, Sec. III.B)")
                echo()

            from optim_wien.converge import run_convergence as _run_conv
            try:
                conv_engine = _run_conv(
                    structure=structure,
                    rmt_result=rmt_result,
                    kmesh_result=kmesh_result,
                    mixing_result=mixing_result,
                    core_valence_result=core_valence_result,
                    gmax_result=gmax_result,
                    lmax_result=lmax_result,
                    basename=basename,
                    work_dir=os.path.abspath(args.output),
                    converge_params=converge_params_set,
                    etol_mRy_per_atom=args.etol,
                    for_forces=for_forces,
                    cluster_submit=args.cluster_submit,
                    quiet=quiet,
                )

                # Override table values with converged values
                if "rkmax" in converge_params_set and conv_engine.final_rkmax > 0:
                    rkmax_result.rkmax = conv_engine.final_rkmax
                if "kmesh" in converge_params_set and conv_engine.final_kmesh[0] > 0:
                    kmesh_result.mesh = conv_engine.final_kmesh
                    kmesh_result.total_points = (
                        conv_engine.final_kmesh[0] *
                        conv_engine.final_kmesh[1] *
                        conv_engine.final_kmesh[2]
                    )
                if "gmax" in converge_params_set:
                    gmax_result.gmax = conv_engine.final_gmax

                # Generate convergence report
                converge_rpt_path = args.converge_report
                if converge_rpt_path is None:
                    converge_rpt_path = os.path.join(
                        args.output, f"{basename}_convergence_report.md"
                    )
                md_report = build_convergence_report(
                    conv_engine, structure, struct_path=args.struct_file,
                    calc_type=args.calc_type,
                )
                with open(converge_rpt_path, "w") as f:
                    f.write(md_report)

                if not quiet:
                    echo()
                    if conv_engine.converged:
                        success(f"Converged. Report: {converge_rpt_path}")
                    else:
                        warn(f"Not fully converged — see report: {converge_rpt_path}")
                    for w in conv_engine.warnings:
                        warn(w)
                    for e in conv_engine.errors:
                        warn(f"  ERROR: {e}")

            except RuntimeError as e:
                if not quiet:
                    warn(f"Convergence engine failed: {e}")
                conv_engine = ConvergenceEngineResult(
                    converged=False,
                    errors=[str(e)],
                )

    elif args.auto_converge:
        if not wien2k_available():
            if not quiet:
                warn("WIEN2k not found in PATH. Cannot auto-converge.")
                echo(style("       Using recommended parameters only.", fg="bright_black"))
            conv_result = ConvergenceResult(
                warnings=["WIEN2k not available — using recommended parameters."]
            )
        else:
            if not quiet:
                header("AUTO-CONVERGENCE", "bright_magenta")
                info("Method", "Blaha hierarchy: k-mesh → RKMAX")
                info("Threshold", "ΔE < 0.1 mRy")
                echo()

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

            if not quiet:
                header("CONVERGENCE RESULTS", "bright_cyan")
                info("Converged", style("YES" if conv_result.converged else "NO",
                       fg="green" if conv_result.converged else "yellow"))
                if conv_result.final_rmts:
                    rmt_str = ", ".join(
                        f"{structure.atoms[i].element}={conv_result.final_rmts[i]:.4f}"
                        for i in range(len(conv_result.final_rmts))
                    )
                    info("Final RMTs", rmt_str)
                if conv_result.final_kmesh[0] > 0:
                    info("Final k-mesh",
                         f"{conv_result.final_kmesh[0]}×{conv_result.final_kmesh[1]}×{conv_result.final_kmesh[2]}")
                info("Final RKMAX", f"{conv_result.final_rkmax}")
                info("Runtime", f"{conv_result.total_runtime:.1f}s")
                for w in conv_result.warnings:
                    warn(w)
                echo()

    runtime = time.time() - t_start

    if not quiet:
        _show_results(rmt_result, rkmax_result, gmax_result, lmax_result,
                      kmesh_result, mixing_result, core_valence_result,
                      structure, conv_result, basename, args, gen_files, runtime)

    has_critical = any("CRITICAL" in w for w in rmt_result.core_leakage_warnings)
    has_final = any("FINAL WARNING" in w for w in rmt_result.overlap_warnings)

    if interactive and not quiet:
        echo()
        _post_run_menu(args, config, structure, rmt_result, rkmax_result,
                       kmesh_result, core_valence_result, basename, has_critical)

    return _build_return(has_critical, has_final, args.quiet,
                         kmesh_result, rkmax_result, core_valence_result)


def _show_results(rmt_result, rkmax_result, gmax_result, lmax_result,
                  kmesh_result, mixing_result, core_valence_result,
                  structure, conv_result, basename, args, gen_files, runtime):
    header("OPTIMIZATION RESULTS", "bright_green")

    # Quick summary card
    rmt_str = ", ".join(
        f"{structure.atoms[i].element}={rmt_result.rmt_values[i]:.3f}"
        for i in range(len(structure.atoms))
    )

    n1, n2, n3 = kmesh_result.mesh

    card_lines = (
        f"  {style('RMT:', fg='bright_yellow', bold=True)}      {style(rmt_str, fg='white')}\n"
        f"  {style('RKMAX:', fg='bright_yellow', bold=True)}    "
        f"{style(f'{rkmax_result.rkmax}', fg='bright_magenta')}    "
        f"{style('GMAX:', fg='bright_yellow', bold=True)}     "
        f"{style(f'{gmax_result.gmax}', fg='bright_cyan')}\n"
        f"  {style('Ecut:', fg='bright_yellow', bold=True)}     "
        f"{style(f'{core_valence_result.ecut:.1f} Ry', fg='bright_cyan')}   "
        f"{style('Mixing:', fg='bright_yellow', bold=True)}  "
        f"{style(f'{mixing_result.scheme} {mixing_result.mixing_factor:.2f}', fg='bright_cyan')}\n"
        f"  {style('k-mesh:', fg='bright_yellow', bold=True)}   "
        f"{style(f'{n1}×{n2}×{n3} ({kmesh_result.total_points} pts)', fg='bright_magenta')}   "
        f"{style('TEMP:', fg='bright_yellow', bold=True)}     "
        f"{style(f'{mixing_result.temp:.4f} Ry', fg='bright_cyan')}\n"
        f"  {style('System:', fg='bright_yellow', bold=True)}   "
        f"{style(kmesh_result.system_type, fg='white')}          "
        f"{style('Runtime:', fg='bright_yellow', bold=True)}  "
        f"{style(f'{runtime:.1f}s', fg='white')}"
    )
    box(card_lines, title="Summary", fg="bright_green")

    # Output files
    echo()
    section("Output Files")
    info("Directory", args.output)
    info("Report", f"{basename}_optimization_report.txt")
    info("Struct", f"{basename}.struct_optimized")
    if not args.no_input_files:
        for name in gen_files:
            info(name, basename)

    # WIEN2k command
    echo()
    section("Ready for WIEN2k")

    prec_flag = getattr(args, "prec", None)
    prec_str = f" -prec {prec_flag}" if prec_flag is not None else ""
    init_cmd = (f"init_lapw -b -rkmax {rkmax_result.rkmax} "
                f"-numk {n1*n2*n3} "
                f"-ecut {int(abs(core_valence_result.ecut))}{prec_str}")
    echo(style(f"    $ {init_cmd}", fg="bright_black"))
    echo(style(f"    $ run_lapw -p", fg="bright_black"))
    echo()

    if getattr(args, "machines", 0) > 0:
        info("Machines", f"{basename}.machines ({args.machines} processes)")

    if getattr(args, "submit", None) == "slurm":
        info("SLURM", f"submit_{basename}_slurm.sh")
    elif getattr(args, "submit", None) == "pbs":
        info("PBS", f"submit_{basename}_pbs.pbs")

    if getattr(args, "vxc", None) == "hse":
        info("Hybrid", f"{basename}.in2c (HSE functional)")

    echo()


def _post_run_menu(args, config, structure, rmt_result, rkmax_result,
                   kmesh_result, core_valence_result, basename, has_critical):
    """Interactive menu after optimization completes."""
    while True:
        idx = menu(
            "What would you like to do next?",
            [
                {"label": "View Full Report", "desc": "Display the complete scientific report"},
                {"label": "View Generated Files", "desc": "List all output files"},
                {"label": "Re-run with Different Settings", "desc": "Go back to parameter selection"},
                {"label": "Edit RMT Manually", "desc": "Adjust one or more RMT values"},
                {"label": "Edit RKMAX", "desc": "Change the RKMAX value"},
                {"label": "Exit", "desc": "Done — continue to WIEN2k"},
            ],
            prompt="Action",
            default=5,
        )

        if idx < 0 or idx == 5:
            break

        if idx == 0:
            rpt_path = os.path.join(args.output,
                                    f"{basename}_optimization_report.txt")
            if os.path.isfile(rpt_path):
                with open(rpt_path) as f:
                    content = f.read()
                clear()
                header("FULL SCIENTIFIC REPORT", "bright_cyan")
                echo(content)
            else:
                warn("Report file not found.")
            input(style("\n  Press Enter to continue...", fg="bright_black"))

        elif idx == 1:
            clear()
            header("GENERATED FILES", "bright_cyan")
            out = args.output
            if os.path.isdir(out):
                for f in sorted(os.listdir(out)):
                    full = os.path.join(out, f)
                    size = os.path.getsize(full)
                    if f.endswith(".struct_optimized"):
                        kind = style("STRUCT", fg="green")
                    elif f.endswith("_report.txt"):
                        kind = style("REPORT", fg="yellow")
                    elif f.endswith(".in0"):
                        kind = style("IN0   ", fg="cyan")
                    elif f.endswith(".in1"):
                        kind = style("IN1   ", fg="cyan")
                    elif f.endswith(".in2"):
                        kind = style("IN2   ", fg="cyan")
                    elif f.endswith(".inm"):
                        kind = style("INM   ", fg="cyan")
                    elif f.endswith(".klist"):
                        kind = style("KLIST ", fg="cyan")
                    else:
                        kind = style("OTHER ", fg="bright_black")
                    echo(f"  {kind}  {f}  ({size} B)")
            input(style("\n  Press Enter to continue...", fg="bright_black"))

        elif idx == 2:
            wizard = InteractiveWizard(struct_file=args.struct_file)
            new_config = wizard.run()
            if new_config:
                config.update(new_config)
                clear()
                banner()
                _run_optimization(config, interactive=True)
                return

        elif idx == 3:
            clear()
            header("MANUAL RMT EDIT", "bright_yellow")
            echo(style("  Current RMT values:", fg="bright_black"))
            for i, a in enumerate(rmt_result.rmt_values):
                echo(f"    [{i + 1}] {structure.atoms[i].element}: {a:.5f} bohr")
            echo()
            echo(style("  Enter new value or leave blank to keep:", fg="bright_black"))
            for i in range(len(rmt_result.rmt_values)):
                a = rmt_result.rmt_values[i]
                raw = input(
                    f"    {structure.atoms[i].element} [{a:.5f}]: "
                ).strip()
                if raw:
                    try:
                        rmt_result.rmt_values[i] = float(raw)
                    except ValueError:
                        warn(f"Invalid number for {structure.atoms[i].element}")
            echo()
            success("RMT values updated. Re-run to regenerate files with new values.")

        elif idx == 4:
            clear()
            header("MANUAL RKMAX EDIT", "bright_yellow")
            info("Current RKMAX", f"{rkmax_result.rkmax}")
            raw = input(style("\n  New RKMAX: ", fg="bright_black")).strip()
            if raw:
                try:
                    rkmax_result.rkmax = float(raw)
                    success(f"RKMAX set to {rkmax_result.rkmax}")
                except ValueError:
                    warn("Invalid number.")

        clear()


def _build_return(has_critical, has_final, quiet, kmesh_result,
                  rkmax_result, core_valence_result):
    class Ret:
        pass
    r = Ret()
    r.exit_code = 0
    if has_critical or has_final:
        if not quiet:
            warn("Critical warnings — please review RMT values manually.")
        r.exit_code = 2
    return r


def _find_struct_file():
    """Auto-detect case.struct files in the current directory."""
    import glob
    candidates = sorted(glob.glob("*.struct"))
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    echo()
    header("MULTIPLE STRUCT FILES FOUND", "bright_yellow")
    echo()
    for i, f in enumerate(candidates):
        echo(f"  {style(f'[{i + 1}]', fg='bright_magenta')} "
             f"{style(f, fg='white')}")
    echo()
    while True:
        try:
            raw = input(f"  Pick one [1-{len(candidates)}]: ").strip()
            idx = int(raw) - 1
            if 0 <= idx < len(candidates):
                return candidates[idx]
        except (ValueError, EOFError, KeyboardInterrupt):
            pass
        echo(style(f"  Enter 1-{len(candidates)}", fg="red"))


def main():
    args = parse_args()

    if args.interactive or (not args.struct_file and sys.stdin.isatty()):
        run_interactive()
        return
    elif args.interactive:
        run_interactive()
        return

    if not args.struct_file:
        detected = _find_struct_file()
        if detected:
            echo()
            echo(style(f"  {BOX['bullet']}  Using: {detected}", fg="bright_green"))
            args.struct_file = detected
        else:
            error("No struct file found in current directory.\n"
                  "  Provide one: opt_wien2k case.struct\n"
                  "  Or run interactive: opt_wien2k -i")

    config = {
        "struct_file": args.struct_file,
        "calc_type": args.calc_type,
        "precision": args.precision,
        "refinement": args.refinement,
        "system_type": args.system_type,
        "bandgap": args.bandgap,
        "vxc": args.vxc,
        "magnetic": args.magnetic,
        "auto_converge": args.auto_converge,
        "converge": args.converge,
        "etol": args.etol,
        "cluster_submit": args.cluster_submit,
        "converge_report": args.converge_report,
        "output": args.output,
        "no_input_files": args.no_input_files,
        "quiet": args.quiet,
        "only": args.only,
    }
    ret = _run_optimization(config, interactive=False)
    sys.exit(ret.exit_code if ret else 0)


if __name__ == "__main__":
    main()
