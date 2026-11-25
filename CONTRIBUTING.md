# Contributing

Contributions are welcome: bug reports, new tools, support for other Ventsim versions, and
solver improvements.

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

## Reporting a bug

Include the Ventsim version (Help > About), the tool you called and its arguments, and the
error or the wrong output. If you can share a model, strip anything confidential first.
