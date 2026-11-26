# Roadmap

ventsim-mcp reads Ventsim `.vsm` models and has its own experimental solvers for airflow,
heat, gas and escape routes. The next step is less about new tools and more about trust:
checking what exists against Ventsim re-runs, then filling the physics that is still missing.

Every task below is a GitHub issue. Pick one, comment that you are on it, and follow
[CONTRIBUTING.md](CONTRIBUTING.md).

## Principles

1. **Check before you extend.** A result is only "matches Ventsim" when a Ventsim re-run says
   so. Reproducing a saved model is a reading check, not an accuracy figure.
2. **Say what is approximate.** Every experimental tool states in its docstring what it
   models, what it leaves out and how it was checked.
3. **No models in the repository.** Ventsim's examples and client models are not ours to
   share. Tests use small synthetic `.vsm` files; results from real models go into issues as
   numbers.
4. **Small pull requests.** One topic, one test, before/after bench output for solver changes.

## Milestones

### 0.15: validation

Check the tools that exist against Ventsim before adding more physics.

| Issue | Task | Needs Ventsim |
|---|---|---|
| #24 | Fan speed response on Tutorial2 | yes |
| #25 | Rock-heat response (`simulate_heat rock_temp`) | yes |
| #26 | Open and re-simulate every kind of written copy | yes |
| #27 | More before/after pairs for `bench/predict.py` | yes |
| #28 | `.vsm` files from other Ventsim versions | a file from that version |
| #29 | Fan pressure at the duty point (~0.7 % median) | no |
| #44 | Example prompts and outputs for every tool | no |
| #45 | Lint and a synthetic-network benchmark in CI | no |
| #46 | Validation report page | no |

### 0.16: HVAC and transient

| Issue | Task |
|---|---|
| #32 | Solve airflow in HVAC models (rooms and doors) |
| #33 | Leakage along auxiliary ventilation ducts |
| #34 | Time-dependent contaminant transport (gas, smoke, DPM) |

### Backlog

| Issue | Task | Area |
|---|---|---|
| #30 | Couple natural ventilation with the heat solution | solver |
| #31 | Speed up the solver on large models | solver |
| #35 | Transient fire: heat release, throttling, flow reversal (after #34) | heat & gas |
| #36 | Rock heat: wet rock and rock age | heat & gas |
| #37 | `shortest_path`: horizontal, duct-free, follow-airflow modes | tools |
| #38 | `goal_seek`: resistance / regulator, velocity / pressure | tools |
| #39 | `sensitivity`: heat mode | tools |
| #40 | `radon_steady`: progeny and dose units | tools |
| #41 | Import airways from DXF | tools |
| #42 | Compare survey measurements (VentLOG / CSV) with the model | tools |
| #43 | Airway size economics (NPV) | tools |

Not planned: Ventsim's 3D drawing interface, menus and toolbars. An MCP server works on the
model file; drawing stays in Ventsim or CAD (`export_dxf`, `render_network`).

## Labels

| Label | Meaning |
|---|---|
| `good first issue` | Self-contained, no Ventsim licence needed, clear test |
| `help wanted` | We would especially like outside help |
| `needs-ventsim` | Checking the result needs a licensed Ventsim DESIGN |
| `area: solver` | Airflow network solver |
| `area: heat & gas` | Heat, moisture, gas, DPM, radon, fire |
| `area: tools` | MCP tools and analysis features |
| `area: validation` | Comparison against Ventsim results |
| `area: file format` | `.vsm` parsing and writing |
| `area: docs & infra` | Documentation, CI, packaging |

## How work flows

```
issue (agree on the approach)
  -> branch from main
  -> code + test in tests/ (synthetic .vsm)
  -> pytest, and bench/fit.py before/after for solver changes
  -> pull request (template), CI green
  -> review, squash merge
  -> a maintainer bumps the version and creates a GitHub release
     (publish.yml sends it to PyPI and the MCP registry)
```

### Contributing without code

If you have a Ventsim licence you can move the project forward without writing Python:
take a `needs-ventsim` issue, make the copy with the tool it names, re-simulate it in Ventsim,
and post the numbers with the "Validation result" issue form. That is the most useful
contribution right now.
