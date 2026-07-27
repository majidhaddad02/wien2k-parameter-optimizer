#!/usr/bin/env bash
#===============================================================================
# WIEN2k Parameter Optimizer — Installer
#===============================================================================
#
# Usage:
#   ./install.sh                  # install to user ~/.local
#   ./install.sh --user           # same as above
#   ./install.sh --prefix /opt    # install to /opt
#   ./install.sh --here           # install in current directory (no pip needed)
#   ./install.sh --system         # system-wide (requires root)
#   ./install.sh --dry-run        # show what would be done
#
# After install, run from anywhere:
#   opt_wien2k case.struct
#   opt_wien2k -i
#
# Requires: python3
#===============================================================================

set -euo pipefail

RED='\033[31m'; GREEN='\033[32m'; YELLOW='\033[33m'; BLUE='\033[34m'
CYAN='\033[36m'; BOLD='\033[1m'; NC='\033[0m'

ok()   { echo -e "  ${GREEN}✓${NC} $*"; }
warn() { echo -e "  ${YELLOW}⚠${NC} $*"; }
err()  { echo -e "  ${RED}✗${NC} $*" >&2; exit 1; }
info() { echo -e "  ${BLUE}→${NC} $*"; }

MODE="user"
PREFIX=""
DRY_RUN=0
FORCE=0

usage() {
    cat <<EOF
${BOLD}WIEN2k Parameter Optimizer — Installer${NC}

Usage:
  ./install.sh [options]

Options:
  --user                Install for current user only (default)
  --prefix PATH         Install to PATH/bin  (e.g. --prefix /opt/local)
  --system              Install system-wide (requires root / sudo)
  --here                Install alias in current directory (no pip)
  --uninstall           Remove the installation
  --dry-run             Show what would be done without doing it
  --force               Overwrite existing installation
  -h, --help            Show this help

After install, use from anywhere:
  ${CYAN}opt_wien2k case.struct${NC}
  ${CYAN}opt_wien2k -i${NC}

EOF
    exit 0
}

# --- Parse args ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --user)     MODE="user" ;;
        --system)   MODE="system" ;;
        --here)     MODE="here" ;;
        --prefix)   MODE="prefix"; PREFIX="$2"; shift ;;
        --dry-run)  DRY_RUN=1 ;;
        --force)    FORCE=1 ;;
        --uninstall) MODE="uninstall" ;;
        -h|--help)  usage ;;
        *)          err "Unknown option: $1. Use --help." ;;
    esac
    shift
done

# --- Find repo root ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$SCRIPT_DIR"

info "Installer directory: $REPO_DIR"

# --- Check python ---
PYTHON=""
for p in python3 python; do
    hash "$p" 2>/dev/null && { PYTHON="$p"; break; }
done
[[ -z "$PYTHON" ]] && err "python3 not found in PATH."

PY_VER=$("$PYTHON" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
ok "Found Python $PY_VER"

# --- Uninstall ---
if [[ "$MODE" == "uninstall" ]]; then
    info "Uninstalling..."

    TARGETS=()
    # pip uninstall
    if hash pip3 2>/dev/null; then
        pip3 uninstall -y opt-wien2k 2>/dev/null && ok "Removed pip package" || true
    fi

    # Shell aliases
    for f in "$HOME/.local/bin/opt_wien2k" "$HOME/.config/opencode/opt_wien2k" \
             "/usr/local/bin/opt_wien2k" "$REPO_DIR/opt_wien2k"; do
        [[ -f "$f" || -L "$f" ]] && rm -f "$f" && ok "Removed $f" || true
    done

    ok "Uninstall complete."
    exit 0
fi

# --- Dependencies check ---
DEPS_OK=1
for mod in setuptools pip; do
    "$PYTHON" -c "import $mod" 2>/dev/null || {
        warn "$mod not installed. Will try without it."
        DEPS_OK=0
    }
done

# --- Mode: --here (simple alias in current dir) ---
if [[ "$MODE" == "here" ]]; then
    ALIAS_FILE="$REPO_DIR/opt_wien2k"
    if [[ -f "$ALIAS_FILE" && "$FORCE" -ne 1 ]]; then
        warn "$ALIAS_FILE already exists. Use --force to overwrite."
    else
        cat > "$ALIAS_FILE" << 'SHELLEOF'
#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec python3 "$SCRIPT_DIR/optimize_wien2k.py" "$@"
SHELLEOF
        chmod +x "$ALIAS_FILE"
        ok "Created: $ALIAS_FILE"
    fi

    echo ""
    echo -e "  ${BOLD}Activate now:${NC}"
    echo -e "    ${CYAN}export PATH=\"$REPO_DIR:\$PATH\"${NC}"
    echo ""
    echo -e "  ${BOLD}Or add to your shell rc:${NC}"
    echo -e "    ${CYAN}echo 'export PATH=\"$REPO_DIR:\$PATH\"' >> ~/.bashrc${NC}"
    echo ""

    # Test
    if [[ "$DRY_RUN" -ne 1 ]]; then
        export PATH="$REPO_DIR:$PATH"
        [[ -x "$ALIAS_FILE" ]] && ok "Ready: opt_wien2k --help" || true
    fi
    exit 0
fi

# --- Mode: pip install ---
if [[ "$MODE" == "system" ]]; then
    if [[ "$EUID" -ne 0 ]]; then
        pip3 install "$REPO_DIR" 2>/dev/null && ok "System install via pip" || {
            warn "Need root. Try: sudo ./install.sh --system"
            exit 1
        }
    else
        pip3 install "$REPO_DIR" 2>/dev/null && ok "System install via pip" || \
            err "System install failed."
    fi
elif [[ "$MODE" == "user" ]]; then
    pip3 install --user "$REPO_DIR" 2>/dev/null && {
        ok "User install via pip"
        info "Binary: ~/.local/bin/opt_wien2k"
        info "Ensure ~/.local/bin is in your PATH"
    } || {
        warn "pip user install failed. Falling back to manual install."

        # Manual: copy to ~/.local/bin
        mkdir -p "$HOME/.local/bin"
        ALIAS_FILE="$HOME/.local/bin/opt_wien2k"
        if [[ -f "$ALIAS_FILE" && "$FORCE" -ne 1 ]]; then
            warn "$ALIAS_FILE already exists. Use --force to overwrite."
        else
            cat > "$ALIAS_FILE" << ALIASEOF
#!/usr/bin/env bash
PYTHONPATH="\$PYTHONPATH:$REPO_DIR"
exec python3 "$REPO_DIR/optimize_wien2k.py" "\$@"
ALIASEOF
            chmod +x "$ALIAS_FILE"
            ok "Manual install: $ALIAS_FILE"
        fi
    }
elif [[ "$MODE" == "prefix" ]]; then
    pip3 install --prefix="$PREFIX" "$REPO_DIR" 2>/dev/null && \
        ok "Install to $PREFIX via pip" || \
        err "pip install --prefix failed."
fi

echo ""
echo -e "  ${BOLD}${GREEN}Install complete.${NC}"
echo ""
echo -e "  ${BOLD}Usage:${NC}"
echo -e "    ${CYAN}opt_wien2k case.struct${NC}           # optimize all parameters"
echo -e "    ${CYAN}opt_wien2k -i${NC}                    # interactive wizard"
echo -e "    ${CYAN}opt_wien2k case.struct --only rmt${NC} # single step"
echo -e "    ${CYAN}opt_wien2k --help${NC}                # see all options"
echo ""
