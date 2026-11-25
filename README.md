# ventsim-mcp

<!-- mcp-name: io.github.agenticeram/ventsim-mcp -->

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
pip install ventsim-mcp
# or the latest code: pip install git+https://github.com/agenticeram/ventsim-mcp
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

`what_if` solves node pressures for the airways of the simulated stage, balancing mass
flow (air expands as it warms) with `dP = R*Q|Q| - Ps(Q) - NVP` on every airway. What it
reads from the file:

- **Stage**: `Stage` is a bit mask and an airway that differs between stages is stored once
  per stage, on top of each other. The stage whose saved airflow balances is used.
- **Junctions**: airway ends closer than 5 cm are one node.
- **Resistance**: the effective value from the saved pressure drop (includes fittings,
  dampers, orifices); scaling in `what_if` multiplies it. Closed airways carry no air.
- **Fans**: density-corrected curves, static or total pressure as the saved results show,
  speed by the affinity laws, curve calibrated through the saved duty point. A fan whose
  saved pressure drop is resistance only is treated as stopped.
- **Boundaries**: portals and dead ends are held at the pressure that reproduces the saved
  flow.
- **Unitemised pressure terms**: natural ventilation pressure, velocity-pressure changes and
  Ventsim's own convergence residual per airway are frozen at their saved values.

Because these inputs come from the saved state, re-solving an unchanged model reproduces
Ventsim's airflow almost exactly. That checks the network was read correctly (stages,
junctions, closed airways, fans); it is not an accuracy figure for changes
(total absolute flow error / total flow, `python bench/fit.py <folder>`):

| Ventsim DESIGN 6 example | Airways solved | Error |
|---|---|---|
| diamondmine | 1,057 | 0.06 % |
| Tutorial2 | 184 | 0.04 % |
| InclinedHeatMine | 810 | 0.02 % |
| the other 17 examples (CoalMine 3,403 airways, Blue_sky, Fire, HVAC, ...) | | < 0.01 % |

**Accuracy for changes.** `Fire` and `Fire2` are the same network with two airway
resistances changed. Predicting Fire2 from Fire plus those two changes
(`python bench/predict.py Fire.vsm Fire2.vsm`) gives an airflow error of 0.03 % of total
flow; it reproduces 96.7 % of the airflow change Ventsim computed. This is the only such
pair among the shipped examples, so treat it as one data point.

Heat, density and fire effects are not re-evaluated after a change. Treat results as an
estimate and confirm in Ventsim.

## Limits

- Read-only on results: simulation, heat and fire modelling stay in Ventsim.
- Editing is limited to text fields; geometry edits would leave results stale.
- Tested on the 20 example models shipped with Ventsim DESIGN 6 (up to 6,700 airways). Other versions untested.

## Development

```
pip install -e ".[test]"
pytest
```

MIT licence.
