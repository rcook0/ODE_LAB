# ODE Lab v3

Includes:
- Python models + integrators + registry + presets: `python/ode_models.py`
- Tiny CLI: `python/ode_lab_cli.py`
- Julia mirror: `julia/ode_models.jl`
- Rust/WASM bindings: `rust-wasm/`
- Web UI: `web/` (live params, export, scenes, 3D attractor projection)

## Build & run UI

```bash
cd rust-wasm
wasm-pack build --target web --release
cd ..
cp -r rust-wasm/pkg web/
python -m http.server -d web 8080
# open http://localhost:8080
```

## Run CLI

```bash
cd python
python -m ode_lab_cli list-models
python -m ode_lab_cli list-presets --preset-file presets.json
python -m ode_lab_cli run --preset duffing_chaotic --preset-file presets.json --out-csv out.csv
```

## v4 additions
- True vector SVG export (SVG primitives)
- Trajectory/attractor caching (view/density changes don't re-integrate)
- Adaptive RK45 for Lorenz/Rössler in WASM
- Advanced systems: Chua + Thomas in WASM; Chen/Halvorsen/Hindmarsh–Rose/Robertson stiff in Python

## v5 step1 (integrator refactor)
- Added `integrate_model(model, params, y0, t0, t1, dt)` (Vec-based RK4) as a generic N-D entry point.
- Refactored existing WASM exports to call the generic integrator internally (UI remains compatible).

## v5 step2 (standard series shape)
- Standardized solver outputs to `{t, y(flat), dim}` across WASM + Python.
- Added `integrate_series()` + `to_series()` in Python.
- UI normalizes any legacy nested `y` arrays defensively.

## Dev helper
- WSL2/Ubuntu: `./scripts/dev.sh setup` then `./scripts/dev.sh all`
- Docs: `docs/dev_wsl2.md`

## v5 Step 3 (Projection pipeline)
- Added projection layer to render x_i vs x_j for any N-D series.
- Default projection (0,1); toggle with keyboard (p) and adjust indices.

## v5 Step 4 (nD registry + metadata)
- Added MODEL_REGISTRY describing dim, params, defaults, and 2D vector field slices.
- Model dropdown and parameter UI are now generated from registry.
- Projection controls added (choose i/j indices).
- ODE integration path prefers `integrate_model()` generic entry point.
