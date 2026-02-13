#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WEB_DIR="$ROOT_DIR/web"
RUST_DIR="$ROOT_DIR/rust-wasm"
PY_DIR="$ROOT_DIR/python"

usage() {
  cat <<'EOF'
Usage: ./scripts/dev.sh [command]

Commands:
  setup       Install common deps in WSL2/Ubuntu (rustup, wasm-pack, venv) - safe to re-run
  wasm        Build wasm (release) and sync to web/pkg
  web         Serve the web UI at http://localhost:8080
  py-venv     Create Python venv + install requirements
  cli [args]  Run the Python CLI inside the venv (auto-creates if missing)
  all         wasm + web (most common loop)

Examples:
  ./scripts/dev.sh setup
  ./scripts/dev.sh all
  ./scripts/dev.sh cli list
  ./scripts/dev.sh cli run --preset lorenz_classic --out /tmp/lorenz.json
EOF
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "Missing command: $1"; exit 1; }
}

setup() {
  if command -v apt >/dev/null 2>&1; then
    echo "[setup] Installing Ubuntu packages..."
    sudo apt update
    sudo apt install -y build-essential curl pkg-config python3 python3-venv python3-pip nodejs npm
  else
    echo "[setup] apt not found; skipping OS packages."
  fi

  if ! command -v cargo >/dev/null 2>&1; then
    echo "[setup] Installing rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
  fi
  # shellcheck disable=SC1091
  source "$HOME/.cargo/env" || true

  if ! command -v wasm-pack >/dev/null 2>&1; then
    echo "[setup] Installing wasm-pack..."
    cargo install wasm-pack
  fi

  py_venv
  echo "[setup] Done."
}

wasm() {
  require_cmd wasm-pack
  echo "[wasm] Building..."
  cd "$RUST_DIR"
  wasm-pack build --target web --release
  cd "$ROOT_DIR"
  echo "[wasm] Sync web/pkg"
  rm -rf "$WEB_DIR/pkg"
  cp -r "$RUST_DIR/pkg" "$WEB_DIR/pkg"
  echo "[wasm] OK"
}

web() {
  require_cmd python3
  echo "[web] Serving $WEB_DIR on http://localhost:8080"
  python3 -m http.server 8080 -d "$WEB_DIR"
}

py_venv() {
  require_cmd python3
  if [ ! -d "$PY_DIR/.venv" ]; then
    echo "[py] Creating venv..."
    python3 -m venv "$PY_DIR/.venv"
  fi
  # shellcheck disable=SC1091
  source "$PY_DIR/.venv/bin/activate"
  python -m pip install -U pip
  if [ -f "$PY_DIR/requirements.txt" ]; then
    pip install -r "$PY_DIR/requirements.txt"
  else
    echo "[py] requirements.txt not found; skipping."
  fi
  echo "[py] Venv ready."
}

cli() {
  py_venv >/dev/null
  # shellcheck disable=SC1091
  source "$PY_DIR/.venv/bin/activate"
  python "$PY_DIR/ode_lab_cli.py" "$@"
}

all() {
  wasm
  web
}

cmd="${1:-}"
shift || true

case "$cmd" in
  setup) setup ;;
  wasm) wasm ;;
  web) web ;;
  py-venv) py_venv ;;
  cli) cli "$@" ;;
  all) all ;;
  ""|-h|--help) usage ;;
  *) echo "Unknown command: ${cmd:-<none>}"; echo; usage; exit 1 ;;
esac
