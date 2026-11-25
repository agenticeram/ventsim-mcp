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
(`dP = R*Q|Q| - Ps(Q)`; the stored resistance and pressures match this relation exactly).
It reports changes relative to an unmodified solve of the same model, and returns a
`baseline_fit` so you can see how well it reproduces Ventsim's own saved flows. On
DemoMine the mean error is 4.4 % (max 7.9 m3/s). It does not model fixed-flow airways,
fires, heat or natural ventilation pressure: on the Tutorial1 and Fire examples the
baseline is off by 20-90 %, and the tool tells you so. Treat results as a first estimate.

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
