# WSL2 dev quickstart

Recommended: keep the repo inside the WSL filesystem (e.g. `~/code/ode_lab`) for speed.

## One-time
```bash
./scripts/dev.sh setup
```

## Daily loop (build wasm + serve UI)
```bash
./scripts/dev.sh all
```

## Run Python CLI
```bash
./scripts/dev.sh cli list
./scripts/dev.sh cli run --preset lorenz_classic --out /tmp/lorenz.json
```

If localhost forwarding isn't working, use:
```bash
hostname -I | awk '{print $1}'
```
then open `http://<ip>:8080`.
