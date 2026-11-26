# ventsim-mcp

<!-- mcp-name: io.github.agenticeram/ventsim-mcp -->

[![agenticeram/ventsim-mcp MCP server](https://glama.ai/mcp/servers/agenticeram/ventsim-mcp/badges/score.svg)](https://glama.ai/mcp/servers/agenticeram/ventsim-mcp)

An [MCP](https://modelcontextprotocol.io) server that lets AI assistants read and inspect
**Ventsim** mine ventilation models (`.vsm` files).

It works on the model file itself, so Ventsim does not need to be running. It reads the
geometry, fans and the results Ventsim saved in the file, and has its own experimental
airflow, temperature and gas solvers for what-if questions; confirm important results in
Ventsim.

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
| `fan_duty_check` | Where each fan runs on its curve: airflow per fan, stall margin, flags fans in or near the stall region or off the curve |
| `fan_summary` | Fans by model with total `Fan Power` as stored in the file (kW) |
| `compare_models` | Diff two saved versions: added/removed airways, airflow and pressure changes |
| `check_flow_balance` | Airflow mass balance at every junction (finds un-simulated or broken networks) |
| `what_if` | Experimental: re-solve the network after scaling airway resistance, switching fans off or changing fan speed, report airflow changes |
| `size_regulator` | Experimental: resistance, pressure drop and approximate opening area of a regulator that brings one airway to a target airflow |
| `trace_airflow` | Main airflow path through an airway from intake to exhaust, with length and frictional pressure along it |
| `simulate_airflow` | Experimental: solve the airflow from the model inputs alone (no saved results needed), with optional changes |
| `simulate_heat` | Experimental: air temperature through the network after airflow changes or with extra heat sources (diesel, fire, coolers) |
| `contaminant_spread` | Experimental: steady-state gas / smoke concentration downstream of release points |
| `render_network` | Draw the network to SVG (plan or section), coloured by any field, fans marked, paths highlighted |
| `export_csv` | Airways to CSV, any columns |
| `edit_airways` | Experimental: write a **copy** with changed Notes / Group / Branch Name, fan speed or fan on/off, airway size, length or friction factor (resistance follows) |
| `add_airways` | Experimental: write a **copy** with new airways (copied from an existing one, joined to the network where the ends meet) |

`list_airways` takes `units="imperial"` (cfm, fpm, ft, in. w.g., hp).

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

## Simulation from inputs (experimental)

`simulate_airflow` does not use any saved result: resistance (`Resistance Total`), fan
curves at their speed and count, fixed-flow airways, fixed-pressure airways (only when
Ventsim's flag for it is set), closed airways, ends open to surface and dead ends (read from
Ventsim's branch flags). Air is taken as incompressible and natural ventilation is ignored,
unless `use_saved_heat=True` borrows both from the saved results. Compared with the airflow
Ventsim saved (`python bench/simulate.py <folder>`):

| Ventsim DESIGN 6 examples | Error, cold | With saved heat |
|---|---|---|
| CoalMine2, Heat Strategies, Tutorial1, RealMine, DemoMine | 0.0 - 0.2 % | same |
| Tutorial2, DevSolids, RealMine 2 | 0.2 - 1.2 % | same |
| CoalMine, diamondmine, without_ducting, Natural Ventilation Effects | 2 - 4 % | 1.3 - 3.5 % |
| Fire, Fire2, InclinedHeatMine (heat-driven) | 6 - 28 % | 0.6 - 2 % |
| Blue_sky | 8.5 % | 5.4 % |
| HVACDemo, A, B (building ducts with room outlets) | not supported | |

The remaining differences come from what the file does not spell out: heat and natural
ventilation, how Ventsim applies total-pressure fan curves, leakage on ducts.

## Heat and gas (experimental)

`simulate_heat` carries the dry-bulb temperature along the airflow, fully mixed at
junctions. Each airway keeps the sensible heat Ventsim saved for it (rock, diesel, fans,
other sources); the auto-compression part of it is split off and follows the new mass flow.
Extra sources (kW, negative for cooling) can be added; a stopped fan's power is taken out
and a speed change scales it with speed cubed. Changes are reported against the unchanged
model, and `baseline_fit` says how well that matches Ventsim.

Re-running an unchanged model with its saved airflow (`python bench/heat.py <folder>`) is a
reproduction check, not an accuracy figure for changes, because the heat comes from the
same saved results. Outlet dry-bulb temperature, airways above 0.5 kg/s, mass-flow
weighted mean / 95th percentile / worst single airway of |error|:

| Ventsim DESIGN 6 example | Airways | Mean | p95 | Max |
|---|---|---|---|---|
| Heat Strategies | 30 | 0.03 K | 0.05 K | 0.05 K |
| diamondmine | 761 | 0.17 K | 0.31 K | 2.9 K |
| Natural Ventilation Effects | 186 | 0.18 K | 0.47 K | 0.6 K |
| InclinedHeatMine | 731 | 0.18 K | 0.65 K | 14 K |
| Blue_sky | 506 | 0.33 K | 0.55 K | 6.4 K |
| CoalMine | 1,989 | 0.36 K | 1.6 K | 8.8 K |
| Fire, Fire2 (fire airways) | 147 | 1.1 K | 2.7 K | 11.5 K |

There is no independent check yet: Fire2's saved temperatures are the same as Fire's, so
the pair used for airflow says nothing about heat. Moisture (latent heat, wet-bulb) and the
change of rock heat with air temperature are not modelled; treat temperatures after large
changes as an estimate.

`contaminant_spread` gives the steady-state concentration of a gas released in airways
(m3/s of pure gas): `c_out = c_in + G/Q`, mixed by airflow at junctions, recirculation
included. It is checked against exact hand calculations in the tests only: the mine
examples hold time-dependent contaminant results, and the one steady result (dust in
HVACDemoB) is in an HVAC model, which is not supported. Reading the model's own source
inputs (`include_file_sources=True`, `Base Injection` x `Base Concentration`) is an
unverified interpretation of Ventsim's fields and is off by default.

## Editing geometry (experimental)

`edit_airways` accepts `Width`, `Height`, `Length` and `kFactor`. Area and perimeter follow
the airway shape (rectangle, round; other profiles are scaled) and `Resistance Total` is
scaled with the Atkinson relation `R ~ k * P * L / A^3`, so shock and density factors Ventsim
folded into it are kept in proportion. `add_airways` copies an existing airway to new end
points; it does not copy fans, fixed flows or pressures, heat or gas sources. A new end
joins the network where it meets an existing airway end (not at a surface connection,
which it would close); an end that meets nothing is a dead end unless listed in
`open_ends`. Run `simulate_airflow` on the copy for the new airflow. Edited and added copies
have not yet been opened in Ventsim itself; check them there before relying on them.

## Limits

- Fire dynamics, time-dependent contaminant spread, moisture and psychrometrics stay in Ventsim.
- Fan curve pressure at the duty point differs from Ventsim's saved value by about 0.7 % (median); Ventsim applies a correction the file does not hold. Fan airflow is exact.
- Building HVAC examples (HVACDemo, A, B) are not supported: their ducts connect through room objects the airway table does not hold.
- Tested on the 20 example models shipped with Ventsim DESIGN 6 (up to 6,700 airways). Other versions untested.

## Contributing

Issues and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). Questions and ideas go to
[Discussions](https://github.com/agenticeram/ventsim-mcp/discussions).

## Development

```
pip install -e ".[test]"
pytest
```

MIT licence.
