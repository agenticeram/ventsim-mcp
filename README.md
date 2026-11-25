# ventsim-mcp

An [MCP](https://modelcontextprotocol.io) server that lets AI assistants read and inspect
**Ventsim** mine ventilation models (`.vsm` files).

It works on the model file itself, so Ventsim does not need to be running. It does **not**
run the ventilation simulation; it reads the geometry, fans and the results Ventsim saved in
the file.

> Unofficial. Not affiliated with or endorsed by Howden / Chasm Consulting. "Ventsim" is
> their trademark. The `.vsm` format is undocumented; it was read from files saved by
> Ventsim DESIGN 6 and may change between versions.

## Tools

| Tool | What it does |
|---|---|
| `model_summary` | Title, Ventsim version, airway and fan counts, total length, depth range |
| `list_airways` | Filter by name/notes or minimum velocity, sort by any field (airflow, length, pressure...) |
| `get_airway` | Every column of one airway |
| `list_fans` | Fan airways with fan model name, duty airflow, fan pressure, power |
| `velocity_check` | Airways above a maximum (or below a minimum) air velocity |
| `fan_summary` | Fans by model with total `Fan Power` as stored in the file (kW) |
| `compare_models` | Diff two saved versions: added/removed airways, airflow and pressure changes |
| `check_flow_balance` | Airflow mass balance at every junction (finds un-simulated or broken networks) |
| `what_if` | Experimental: re-solve the network after scaling airway resistance or switching fans off, report airflow changes |
| `export_csv` | Airways to CSV, any columns |
| `edit_airways` | Experimental: write a **copy** with changed Notes / Group / Branch Name |

All values are in the SI units Ventsim stores internally (m, m3/s, Pa, kW).

## Install

```
pip install git+https://github.com/agenticeram/ventsim-mcp
```

Claude Desktop / Claude Code config:

```json
{ "mcpServers": { "ventsim": { "command": "ventsim-mcp" } } }
```

Then ask, for example: *"Open C:\models\mine.vsm and list the 10 airways with the highest velocity."*

## File format notes

A `.vsm` is a gzip-compressed, CRLF-separated, tab-delimited text file. The airway table starts
at a header line beginning `MAIN<TAB>Unique No`, which names all ~288 columns; the parser uses
those names, never fixed positions. Reading and writing is byte-exact when nothing is edited. A `.txt` export of the same file (uncompressed) is read too.
The CSV/NBC/XLSX exports are headerless or derived dumps of the same tables and are not supported; use the `.vsm`/`.txt`.

## What-if solver (experimental)

`what_if` solves node pressures with the resistances and fan curves stored in the file
(`dP = R*Q|Q| - Ps(Q) - NVP`). It honours fan speed (affinity laws), fixed-flow airways,
fixed-pressure airways and the natural ventilation pressure saved per airway, and reports
changes relative to an unmodified solve of the same model, together with a `baseline_fit`
showing how well the solver reproduces the flows Ventsim saved. It uses sparse matrices and
solves a 1,300-airway model in about two seconds.

Fan curves are density-corrected, may be static or total pressure (decided from the
saved results of each fan airway), are forced monotone and are extended past their last
point. Flow error against the saved results (total absolute flow error / total flow;
`python bench/fit.py <folder>` reproduces this for any folder of models):

| Example model | Error |
|---|---|
| Tutorial1, RealMine, RealMine 2, Heat Strategies | < 0.1 % |
| DemoMine | 0.13 % |
| without_ducting | 0.6 % |
| Fire, Fire2 | 0.9 % |
| Natural Ventilation Effects | 1.5 % |
| InclinedHeatMine | 1.8 % |
| Blue_sky | 28 %: not reliable (see below) |

Blue_sky links parts of the network that are not joined at a shared point (ducting or
leakage connections); the airflow saved by Ventsim does not balance at 49 junctions when
airways are matched by their end coordinates, so the solver cannot reproduce it. `what_if`
reports `junctions_unbalanced_in_saved_flows` and `reliable` so such models are flagged
instead of returning misleading numbers.

Heat, density and fire effects are not re-evaluated after a change: natural ventilation
pressure stays at the saved value. Treat results as a first estimate and confirm in Ventsim.

## Limits

- Read-only on results: simulation, heat and fire modelling stay in Ventsim.
- Editing is limited to text fields; geometry edits would leave results stale.
- Tested on the example models shipped with Ventsim DESIGN 6 (up to 5,200 airways). Other versions untested.

## Development

```
pip install -e ".[test]"
pytest
```

MIT licence.
