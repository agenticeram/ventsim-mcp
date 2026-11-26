# Contributing

Contributions are welcome: bug reports, new tools, support for other Ventsim versions, and
solver improvements.

## Finding a task

Open tasks are GitHub issues, grouped by milestone in [ROADMAP.md](ROADMAP.md). Look for
`good first issue` (no Ventsim licence needed) or `needs-ventsim` (you can help by
re-simulating in Ventsim, no code). Comment on the issue before you start so two people do
not do the same work.

## Getting started

```
git clone https://github.com/agenticeram/ventsim-mcp
cd ventsim-mcp
pip install -e ".[test]"
pytest
```

## Ground rules

- Open an issue first for larger changes, so we can agree on the approach.
- One topic per pull request. Add or update a test in `tests/` for every change.
- Do not commit Ventsim's example models or any client model: they are not ours to
  redistribute. Tests build small synthetic `.vsm` files instead (see `tests/test_vsm.py`).
- The `.vsm` format is undocumented. When you rely on a column, say in the PR how you
  checked its meaning (which models, what you compared).

## Solver changes

Run the benchmark on the examples that ship with Ventsim (or your own models) before and
after your change and paste both outputs in the PR:

```
python bench/fit.py "C:/Program Files/Ventsim DESIGN 6/examples"
python bench/predict.py Fire.vsm Fire2.vsm
```

`fit.py` must stay under 0.1 % on every example. `predict.py` is the real accuracy check:
if you have two saved versions of the same network that differ by a known change, a new
pair for it is one of the most useful contributions you can make.

## Validation results

If you have a Ventsim licence, the most useful contribution is checking a tool against a
Ventsim re-run. Take a `needs-ventsim` issue, make the copy with the tool it names, open and
re-simulate it in Ventsim, and post the numbers with the "Validation result" issue form.
Post numbers, not the model.

## Pull requests

Branch from `main`, keep the pull request to one topic, fill in the template and make CI
green. Pull requests are squash merged. Versions are bumped and released by a maintainer; a
GitHub release publishes to PyPI and the MCP registry.

## Reporting a bug

Include the Ventsim version (Help > About), the tool you called and its arguments, and the
error or the wrong output. If you can share a model, strip anything confidential first.
