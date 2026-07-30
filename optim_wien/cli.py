"""
Interactive CLI / TUI for WIEN2k Parameter Optimizer.

Provides:
  - Interactive wizard mode (--interactive / -i)
  - Rich colorized output with box-drawing
  - Number-based menus (no external deps)
  - Progress indicators

Usage:
    from optim_wien.cli import Cli, echo
"""

import os
import sys
import time
import shutil

_FG = {
    "black": "\033[30m", "red": "\033[31m", "green": "\033[32m",
    "yellow": "\033[33m", "blue": "\033[34m", "magenta": "\033[35m",
    "cyan": "\033[36m", "white": "\033[37m",
    "bright_black": "\033[90m", "bright_red": "\033[91m",
    "bright_green": "\033[92m", "bright_yellow": "\033[93m",
    "bright_blue": "\033[94m", "bright_magenta": "\033[95m",
    "bright_cyan": "\033[96m", "bright_white": "\033[97m",
}

_BG = {
    "black": "\033[40m", "red": "\033[41m", "green": "\033[42m",
    "yellow": "\033[43m", "blue": "\033[44m", "magenta": "\033[45m",
    "cyan": "\033[46m", "white": "\033[47m",
}

_STYLE = {
    "bold": "\033[1m", "dim": "\033[2m", "italic": "\033[3m",
    "underline": "\033[4m", "blink": "\033[5m", "reverse": "\033[7m",
    "hidden": "\033[8m", "strike": "\033[9m",
}

_RESET = "\033[0m"

# Box-drawing chars
BOX = {
    "tl": "╭", "tr": "╮", "bl": "╰", "br": "╯",
    "h": "─", "v": "│",
    "tl2": "┏", "tr2": "┓", "bl2": "┗", "br2": "┛",
    "h2": "━", "v2": "┃",
    "dot": "●", "arrow": "→", "check": "✓", "cross": "✗",
    "bullet": "•", "dash": "─", "double_dash": "═",
}

_TERM_WIDTH = None


def term_width():
    global _TERM_WIDTH
    if _TERM_WIDTH is None:
        _TERM_WIDTH = shutil.get_terminal_size((80, 24)).columns
    return _TERM_WIDTH


def style(text, fg=None, bg=None, bold=False, dim=False, italic=False,
           underline=False):
    codes = _FG.get(fg, "") if fg else ""
    codes += _BG.get(bg, "") if bg else ""
    if bold:
        codes += _STYLE["bold"]
    if dim:
        codes += _STYLE["dim"]
    if italic:
        codes += _STYLE["italic"]
    if underline:
        codes += _STYLE["underline"]
    if not codes:
        return text
    return f"{codes}{text}{_RESET}"


def echo(text="", **kwargs):
    print(text, **kwargs)


def box(text, title=None, fg="cyan", width=None):
    w = width or min(term_width() - 4, 80)
    inner = w - 2
    echo(style(f"{BOX['tl']}{BOX['h'] * inner}{BOX['tr']}", fg=fg))
    if title:
        title_line = f" {title[:inner-2]} "
        pad = inner - len(title) - len(title_line)
        echo(style(f"{BOX['v']}", fg=fg) + style(title_line, fg=fg, bold=True)
              + style(f"{BOX['v']}", fg=fg))
        echo(style(f"{BOX['v']}{BOX['h'] * inner}{BOX['v']}", fg=fg))
    for line in text.split("\n"):
        visible = strip_ansi(line)
        if len(visible) <= inner:
            line = line + " " * (inner - len(visible))
        echo(style(f"{BOX['v']} ", fg=fg) + line
             + style(f" {BOX['v']}", fg=fg))
    echo(style(f"{BOX['bl']}{BOX['h'] * inner}{BOX['br']}", fg=fg))


def strip_ansi(s):
    import re
    return re.sub(r"\033\[[0-9;]*m", "", s)


def progress_bar(current, total, width=40, label="", fg="green"):
    pct = current / total if total > 0 else 1
    filled = int(width * pct)
    bar = style("█" * filled, fg=fg) + style("░" * (width - filled), fg="bright_black")
    sys.stdout.write(f"\r  {label} {bar} {pct*100:.0f}%")
    sys.stdout.flush()
    if current >= total:
        sys.stdout.write("\n")


def spinner(frames, interval=0.1, message=""):
    """Generator that yields spinner frames."""
    import itertools
    for i, frame in enumerate(itertools.cycle(frames)):
        yield frame


def banner():
    """Display ASCII banner."""
    text = r"""
   ^    ____    ____    __ __   _____   _  _   ___  __ __ 
  <^>  / _  \ /  _  \ |  ^  \ |  _  | | |/ / / _ \|  \\  \
  <^> |  | | ||  |  | ||  ^   || | | | |   / |  __/|      |
  <^> |  |_| /\__/|  ||  |\  || |_| | |  <   \___/|  /\  |
   V   \____/    \_/ |_| | \_|\_____| |_|\_\      |_|  \_|

   """  # noqa: W605
    for line in text.split("\n"):
        echo(style(line, fg="cyan"))


def clear():
    os.system("clear" if os.name != "nt" else "cls")


def header(text, fg="bright_cyan"):
    w = term_width()
    echo()
    echo(style(f"  {text}  ".center(w, BOX["h2"]), fg=fg, bold=True))
    echo()


def section(text, fg="yellow"):
    echo(style(f"  {BOX['dot']} {text}", fg=fg, bold=True))


def info(key, value, indent=2):
    key_str = style(f"{key}:", fg="bright_black", bold=True)
    val_str = style(str(value), fg="white")
    echo(f"{' ' * indent}{key_str} {val_str}")


def warn(text):
    echo(style(f"  ⚠  {text}", fg="yellow"))


def error(text):
    echo(style(f"  ✗  {text}", fg="red", bold=True))
    sys.exit(1)


def success(text):
    echo(style(f"  {BOX['check']}  {text}", fg="green", bold=True))


def menu(title, options, prompt="Choice", default=None, clear_first=False):
    """Display a numbered menu and return the selected index (0-based)."""
    if clear_first:
        clear()
    echo()
    box(title, fg="bright_cyan")
    echo()
    for i, opt in enumerate(options):
        label = opt["label"]
        desc = opt.get("desc", "")
        key = style(f"  [{i + 1}]", fg="bright_magenta")
        lbl = style(label, fg="white", bold=True)
        dsc = style(f" — {desc}", fg="bright_black") if desc else ""
        echo(f"{key} {lbl}{dsc}")
    if default is not None:
        echo(style(f"\n  [{len(options) + 1}] back / cancel", fg="bright_black"))
    echo()

    while True:
        try:
            raw = input(f"  {prompt} "
                        f"{style(f'[1-{len(options)}]', fg='bright_black')}: ").strip()
            if not raw and default is not None:
                return default
            idx = int(raw) - 1
            if 0 <= idx < len(options):
                return idx
            if default is not None and idx == len(options):
                return -1
        except (ValueError, EOFError, KeyboardInterrupt):
            if default is not None:
                return -1
        echo(style(f"  Enter 1-{len(options)}", fg="red"))


def confirm(text, default=True):
    yes = "Y" if default else "y"
    no = "n" if default else "N"
    while True:
        try:
            raw = input(
                f"  {text} [{yes}/{no}]: ").strip().lower()
            if not raw:
                return default
            if raw in ("y", "yes"):
                return True
            if raw in ("n", "no"):
                return False
        except (EOFError, KeyboardInterrupt):
            return default
        echo(style("  Enter y/n", fg="red"))


def struct_input(default_path=None):
    """Prompt for struct file path with validation."""
    clear()
    banner()
    header("STEP 1: Select Structure File")

    while True:
        if default_path:
            hint = f" [default: {default_path}]"
        else:
            hint = ""
        raw = input(f"  Path to case.struct file{hint}: ").strip()
        if not raw and default_path:
            raw = default_path
        if not raw:
            echo(style("  ✗  Please enter a file path", fg="red"))
            continue
        if os.path.isfile(raw):
            return raw
        echo(style(f"  ✗  '{raw}' not found", fg="red"))


class InteractiveWizard:
    """Full interactive wizard for parameter selection."""

    def __init__(self, struct_file=None):
        self.struct_file = struct_file
        self.calc_type = "scf"
        self.precision = "medium"
        self.refinement = "medium"
        self.system_type = None
        self.vxc = "pbe"
        self.magnetic = False
        self.auto_converge = False
        self.output = "./optim_results"
        self.no_input_files = False
        self.quiet = False

    def run(self):
        clear()
        banner()
        echo(style("  WIEN2k Interactive Parameter Optimizer", fg="bright_cyan", bold=True))
        echo(style("  Scientific references: Blaha et al., J. Chem. Phys. 152, 074101 (2020)", fg="bright_black"))
        time.sleep(0.8)

        self.struct_file = self._step_struct()
        if self.struct_file is None:
            return None

        self.calc_type = self._step_calc_type()
        if self.calc_type is None:
            return None

        self.precision = self._step_precision()
        if self.precision is None:
            return None

        self.refinement = self._step_refinement()
        if self.refinement is None:
            return None

        self.vxc = self._step_vxc()
        if self.vxc is None:
            return None

        self.magnetic = self._step_magnetic()
        if self.magnetic is None:
            return None

        self.auto_converge = self._step_auto_converge()
        if self.auto_converge is None:
            return None

        self.output = self._step_output()
        if self.output is None:
            return None

        return self._summary()

    def _step_struct(self):
        clear()
        banner()
        header("STEP 1: Select Structure File")
        echo(style("  Provide the path to your WIEN2k case.struct file.", fg="bright_black"))
        echo(style("  Example: BaTiO3.struct  or  /path/to/case.struct", fg="bright_black"))
        echo()
        while True:
            default = self.struct_file or ""
            hint = f" [default: {default}]" if default else ""
            raw = input(f"  Struct file{hint}: ").strip()
            if not raw and default:
                raw = default
            if raw.lower() in ("q", "quit", "exit", "back"):
                return None
            if not raw:
                echo(style("  ✗  Please enter a file path", fg="red"))
                continue
            if os.path.isfile(raw):
                echo()
                success(f"Found: {raw}")
                time.sleep(0.5)
                return raw
            echo(style(f"  ✗  '{raw}' not found. Try again or type 'back'.", fg="red"))

    def _step_calc_type(self):
        clear()
        banner()
        header("STEP 2: Calculation Type")
        idx = menu(
            "What type of WIEN2k calculation are you running?",
            [
                {"label": "SCF", "desc": "Self-consistent field — energy, DOS, band structure"},
                {"label": "Forces", "desc": "Force convergence and geometry relaxation"},
                {"label": "Relaxation", "desc": "Full volume + internal coordinate relaxation"},
                {"label": "Optimization", "desc": "Optimize internal positions only"},
                {"label": "EOS", "desc": "Equation of state — energy vs volume"},
                {"label": "EFG", "desc": "Electric field gradient calculation"},
            ],
            prompt="Select",
            default=0,
            clear_first=False,
        )
        if idx < 0:
            return None
        return ["scf", "forces", "relaxation", "optimization", "eos", "efg"][idx]

    def _step_precision(self):
        clear()
        banner()
        header("STEP 3: Precision Level")
        idx = menu(
            "How accurate do you need the calculation?",
            [
                {"label": "Screening", "desc": "Quick scan — fast, low accuracy"},
                {"label": "Coarse", "desc": "Rough estimates — minimal cost"},
                {"label": "Medium", "desc": "Production quality — good balance (recommended)"},
                {"label": "High", "desc": "Publication quality — stricter convergence"},
                {"label": "Very High", "desc": "Benchmark quality — maximum accuracy"},
            ],
            prompt="Select",
            default=2,
        )
        if idx < 0:
            return None
        return ["screening", "coarse", "medium", "high", "very_high"][idx]

    def _step_refinement(self):
        clear()
        banner()
        header("STEP 4: k-Mesh Refinement")
        idx = menu(
            "How fine do you want the k-point mesh?",
            [
                {"label": "Coarse", "desc": "Quick scans, ~half the recommended density"},
                {"label": "Medium", "desc": "Standard production (recommended)"},
                {"label": "Fine", "desc": "Double density for accurate energies"},
                {"label": "Very Fine", "desc": "4× density for DOS and optical properties"},
            ],
            prompt="Select",
            default=1,
        )
        if idx < 0:
            return None
        return ["coarse", "medium", "fine", "very_fine"][idx]

    def _step_vxc(self):
        clear()
        banner()
        header("STEP 5: Exchange-Correlation Functional")
        idx = menu(
            "Which XC functional do you use?",
            [
                {"label": "PBE", "desc": "GGA — Perdew-Burke-Ernzerhof (recommended)"},
                {"label": "PBEsol", "desc": "GGA — revised for solids, better lattice constants"},
                {"label": "LDA", "desc": "Local density approximation"},
                {"label": "WC", "desc": "Wu-Cohen GGA"},
                {"label": "SCAN", "desc": "meta-GGA — strongly constrained"},
                {"label": "HSE", "desc": "Heyd-Scuseria-Ernzerhof hybrid functional"},
            ],
            prompt="Select",
            default=0,
        )
        if idx < 0:
            return None
        return ["pbe", "pbesol", "lda", "wc", "scan", "hse"][idx]

    def _step_magnetic(self):
        clear()
        banner()
        header("STEP 6: Magnetic Configuration")
        idx = menu(
            "Is this a spin-polarized (magnetic) calculation?",
            [
                {"label": "No", "desc": "Non-magnetic (most semiconductors, insulators)"},
                {"label": "Yes", "desc": "Spin-polarized — ferromagnetic / antiferromagnetic"},
            ],
            prompt="Select",
            default=0,
        )
        if idx < 0:
            return None
        return idx == 1

    def _step_auto_converge(self):
        clear()
        banner()
        header("STEP 7: Auto-Convergence")
        idx = menu(
            "Enable auto-convergence of k-mesh and RKMAX?",
            [
                {"label": "No", "desc": "Generate recommended parameters only (recommended)"},
                {"label": "Yes", "desc": "Run WIEN2k SCF cycles to converge parameters"},
            ],
            prompt="Select",
            default=0,
        )
        if idx < 0:
            return None
        return idx == 1

    def _step_output(self):
        clear()
        banner()
        header("STEP 8: Output Directory")
        while True:
            raw = input(f"  Output directory [default: {self.output}]: ").strip()
            if not raw:
                return self.output
            if raw.lower() in ("q", "quit", "back"):
                return None
            if os.path.isdir(raw) or not os.path.exists(raw):
                return raw
            echo(style(f"  ✗  '{raw}' exists and is not a directory", fg="red"))

    def _summary(self):
        os.makedirs(self.output, exist_ok=True)
        clear()
        banner()
        header("CONFIGURATION SUMMARY", fg="bright_green")

        box(
            f"  {style('Struct:', fg='bright_black', bold=True)}  "
            f"{style(self.struct_file, fg='white')}\n"
            f"  {style('Calculation:', fg='bright_black', bold=True)}  "
            f"{style(self.calc_type, fg='bright_magenta')}\n"
            f"  {style('Precision:', fg='bright_black', bold=True)}  "
            f"{style(self.precision, fg='bright_yellow')}\n"
            f"  {style('k-Mesh:', fg='bright_black', bold=True)}  "
            f"{style(self.refinement, fg='bright_yellow')}\n"
            f"  {style('XC Functional:', fg='bright_black', bold=True)}  "
            f"{style(self.vxc.upper(), fg='bright_cyan')}\n"
            f"  {style('Magnetic:', fg='bright_black', bold=True)}  "
            f"{style('Yes' if self.magnetic else 'No', fg='bright_cyan')}\n"
            f"  {style('Auto-Converge:', fg='bright_black', bold=True)}  "
            f"{style('Yes' if self.auto_converge else 'No', fg='bright_cyan')}\n"
            f"  {style('Output:', fg='bright_black', bold=True)}  "
            f"{style(self.output, fg='white')}",
            fg="bright_green",
        )

        echo()
        if not confirm("Proceed with these settings?", default=True):
            return None
        return self.to_dict()

    def to_dict(self):
        return {
            "struct_file": self.struct_file,
            "calc_type": self.calc_type,
            "precision": self.precision,
            "refinement": self.refinement,
            "system_type": self.system_type,
            "bandgap": None,
            "vxc": self.vxc,
            "magnetic": self.magnetic,
            "auto_converge": self.auto_converge,
            "converge": None,
            "etol": 0.1,
            "cluster_submit": False,
            "converge_report": None,
            "output": self.output,
            "no_input_files": self.no_input_files,
            "quiet": self.quiet,
            "only": None,
            "machines": 0,
            "submit": None,
        }


class ProgressTracker:
    """Show progress during optimization steps."""

    STEPS = {
        "rmt": ("Optimizing muffin-tin radii", "bright_green"),
        "rkmax": ("Optimizing plane-wave cutoff", "bright_cyan"),
        "gmax": ("Optimizing Fourier expansion", "bright_blue"),
        "lmax": ("Optimizing angular momentum cutoffs", "bright_magenta"),
        "kmesh": ("Generating adaptive k-point mesh", "bright_yellow"),
        "mixing": ("Configuring SCF mixing scheme", "bright_green"),
        "core": ("Setting core-valence separation", "bright_cyan"),
        "report": ("Generating scientific report", "bright_blue"),
        "input_files": ("Generating WIEN2k input files", "bright_magenta"),
    }

    def __init__(self):
        self.current = 0
        self.active = 0

    def start(self, active_steps=None):
        self.current = 0
        self.active = active_steps or len(self.STEPS)
        echo()
        w = term_width()
        echo(style("  OPTIMIZATION PROGRESS  ".center(w, BOX["h2"]),
                   fg="cyan", bold=True))
        echo()

    def step(self, step_name):
        if step_name not in self.STEPS:
            return
        self.current += 1
        desc, fg = self.STEPS[step_name]
        marker = style(BOX["dot"], fg=fg)
        label = style(f"{desc}...", fg=fg)
        count = style(f"[{self.current}/{self.active}]", fg="bright_black")
        sys.stdout.write(f"  {marker} {label} {count}")
        sys.stdout.flush()
        self._last_desc = desc
        self._last_fg = fg

    def done(self, text="OK"):
        desc = getattr(self, "_last_desc", "")
        fg_obj = getattr(self, "_last_fg", "green")
        sys.stdout.write(f"\r{' ' * (term_width() - 2)}")
        sys.stdout.write(f"\r  {style(BOX['check'], fg='green')} "
                         f"{style(desc, fg='green')}"
                         f"  {style(text, fg='green', bold=True)}")
        sys.stdout.write("\n")
        sys.stdout.flush()

    def fail(self, text="FAILED"):
        sys.stdout.write(f"  {style(text, fg='red', bold=True)}\n")
        sys.stdout.flush()

    def finish(self):
        echo()
        echo(style("  All steps complete!  ".center(term_width()),
                   fg="green", bold=True))
        echo()
