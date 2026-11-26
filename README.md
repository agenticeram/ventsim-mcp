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
| `simulate_heat` | Experimental: dry- and wet-bulb temperature and humidity after airflow changes or with extra heat and moisture sources (diesel, fire, coolers) |
| `contaminant_spread` | Experimental: steady-state gas / smoke concentration downstream of release points |
| `escape_routes` | Walking distance from every airway to the surface and to the nearest refuge (Ventsim's own rule), optionally with airways blocked |
| `recirculation_check` | Closed loops in the airflow (fan recirculation and other circulating air) |
| `contaminant_arrival` | Time for smoke or gas to reach each airway from release points |
| `diesel_exhaust` | Experimental: diesel particulate (DPM) concentration and air per kW of diesel |
| `compliance_check` | Airways outside velocity, diesel air, wet-bulb and escape-distance limits (your limits) |
| `fan_select` | Fans from the model's fan library that meet a duty, with speed, efficiency, power and stall margin |
| `energy_report` | Fan energy and yearly cost, and the airways where the air power is lost |
| `calibrate_airflow` | Experimental: fit airway resistances to measured airflows (survey calibration) |
| `list_stages` | Stages in the model; most analysis tools take `stage=` to study another one |
| `write_results` | Experimental: write a **copy** with this server's airflow and temperatures in it |
| `render_network` | Draw the network to SVG (plan or section), coloured by any field, fans marked, paths highlighted |
| `export_csv` | Airways to CSV, any columns |
| `edit_airways` | Experimental: write a **copy** with changed Notes / Group / Branch Name, fan speed or fan on/off, airway size, length or friction factor (resistance follows), parallel airways (Pathways), or an airway closed |
| `add_airways` | Experimental: write a **copy** with new airways (copied from an existing one, joined to the network where the ends meet) |
| `model_check` | Ventsim warnings and errors found from the file: loose ends, duplicates, blocked fans, fans or stoppings in series, conflicting fixed flows |
| `fix_model` | Write a **copy** with the safe fixes: duplicates merged with Pathways, loose ends closed, blocked fans off |

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
resistances changed (x2.26 and x1.41; 25 more differ only in the last digits). Predicting Fire2 from Fire plus those changes
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

Moisture is carried the same way: each airway adds the water Ventsim saved as evaporating
there (latent heat / 2450 kJ/kg), mixed by mass at junctions, and the wet bulb follows from
dry bulb, moisture and pressure (standard psychrometrics; Ventsim's saved latent heat and
wet bulbs agree with it to 0.4 % in the heat examples).

Re-running an unchanged model with its saved airflow (`python bench/heat.py <folder>`) is a
reproduction check, not an accuracy figure for changes, because the heat comes from the
same saved results. Outlet temperature, airways above 0.5 kg/s, mass-flow weighted mean /
95th percentile / worst single airway of |error|:

| Ventsim DESIGN 6 example | Airways | Dry bulb | Wet bulb |
|---|---|---|---|
| Heat Strategies | 30 | 0.03 / 0.05 / 0.05 K | 0.01 / 0.02 / 0.02 K |
| diamondmine | 761 | 0.17 / 0.31 / 2.9 K | 0.04 / 0.11 / 0.9 K |
| Natural Ventilation Effects | 186 | 0.18 / 0.47 / 0.6 K | 0.05 / 0.14 / 0.2 K |
| InclinedHeatMine | 731 | 0.18 / 0.65 / 14 K | 0.08 / 0.80 / 8 K |
| Blue_sky | 506 | 0.33 / 0.55 / 6.4 K | 0.10 / 0.33 / 3.3 K |
| CoalMine | 1,989 | 0.36 / 1.6 / 8.8 K | 0.41 / 1.1 / 7.3 K |
| Fire, Fire2 (fire airways) | 147 | 1.1 / 2.7 / 11.5 K | 3.3 / 7.4 / 11 K |

There is no independent check yet: Fire2's saved temperatures are the same as Fire's, so
the pair used for airflow says nothing about heat. Moisture (latent heat, wet-bulb) and the
change of rock heat with air temperature are not modelled; treat temperatures after large
changes as an estimate.

Heat loads stay as saved: rock heat is not computed from rock properties (Ventsim takes the
virgin rock temperature and rock age from global climate settings for most airways, and
where those live in the file is not known yet).

`contaminant_spread` gives the steady-state concentration of a gas released in airways
(m3/s of pure gas): `c_out = c_in + G/Q`, mixed by airflow at junctions, recirculation
included. It is checked against exact hand calculations in the tests only: the mine
examples hold time-dependent contaminant results, and the one steady result (dust in
HVACDemoB) is in an HVAC model, which is not supported. Reading the model's own source
inputs (`include_file_sources=True`, `Base Injection` x `Base Concentration`) is an
unverified interpretation of Ventsim's fields and is off by default.

`contaminant_arrival` gives the first arrival of smoke or gas carried at the mean air
velocity, the rule behind Ventsim's `Spread Time`: it matches all 131 airways of Tutorial1
to within 0.1 % (`python bench/spread.py <folder>`). InclinedHeatMine's saved spread times
are 1 - 4 % faster than its saved airflow allows, so they look like results of an earlier
run; HVACDemoB is an HVAC model.

`diesel_exhaust` mixes diesel particulate from the engines in the model (kW x an emission
factor you give, g/kWh) along the airflow. Ventsim's own saved DPM values do not follow a steady mixing
model (they look time-dependent), so this tool is not checked against Ventsim.

## Escape routes and recirculation

`escape_routes` reproduces Ventsim's `distancetosurface` and `distancetorefuge` exactly:
shortest walk along airways from the middle of each airway to a surface connection or to the
middle of a refuge airway (`Branch Flags2` bit 5); resistances marked impassable in the
model's resistance table (stoppings, roof falls) cannot be walked. All 3,262 airways of the
two examples that hold these results match (`python bench/escape.py <folder>`).
`recirculation_check` finds the strongly connected parts of the directed airflow. Its fan
loop in Fire2 is exactly the seven airways Ventsim marks as recirculation; it also reports
loops without fans (and large circulating groups in Blue_sky and RealMine) that Ventsim
does not mark, so read those as "air goes round here", not as Ventsim's recirculation.

## Design tools

- `fan_select` reads the fan library in the model (static and total pressure, efficiency,
  power per curve point) and finds the speed at which each fan meets a duty. Checked on the
  running fans of the examples (`python bench/fan_select.py <folder>`): for 99 of 103 fans
  with a library curve (for the other 4 no speed from 20 to 250 % fits), the fan's own curve meets its saved duty, at a speed within
  1.3 points of the saved `Fan Speed` (6 points where only a total-pressure curve is given).
- `energy_report` uses the saved fan power (Ventsim's `Fan Power` = airflow x fan total
  pressure / efficiency) and splits the air power over the airways.
- `compliance_check` limits default to common practice values (0.25 - 6 m/s, 20 m/s in
  shafts, 0.06 m3/s per diesel kW, 28 C wet bulb); they are not any particular regulation,
  so pass the values that apply to you.
- `calibrate_airflow` fits resistance factors to measured airflows by least squares on the
  what_if solver. Using Fire2's airflow in the 27 airways whose stored resistance differs
  from Fire as "measurements" on Fire (`python bench/calibrate.py Fire.vsm Fire2.vsm`)
  brings the whole network to within 0.002 % of Fire2's airflow in 8 s. The resistances
  themselves are not recovered: the two airways that really changed (x2.26 and x1.41) come
  out as x1.52 and x1.08, the rest of the change lands on airways in series with them, which
  the measurements cannot tell apart. Treat the factors as one fit, not as measured values.
- `write_results` writes airflow, end pressures, airway pressure drop, fan pressure and
  temperatures into a copy (what_if solver by default). The other tools read such a copy
  consistently: `what_if` on a written copy reproduces it with 0.00 % error on Fire and
  CoalMine. The changes themselves (resistance factor, Fan Off, Fan Speed) go into the
  inputs too, so re-simulating the copy keeps them. Heat terms and gases are left as they
  were. See "Checked in Ventsim" below for copies re-simulated in Ventsim itself.
- `compliance_check` takes the wet bulb Ventsim saved.
- Stages: `simulate_airflow(stage=...)` solves another stage from its inputs and
  `escape_routes(stage=...)` checks its escape routes; the other tools use the saved
  results, which belong to the saved stage.
- Closing an airway (`edit_airways`, field `Closed`) writes it the way Ventsim stores
  closed airways in the examples: resistance type 1, resistance 1000.

Leakage: in these examples leakage paths are ordinary high-resistance airways, already in
the network solve; nothing else is needed for them.

## Editing geometry (experimental)

`edit_airways` accepts `Width`, `Height`, `Length` and `kFactor`. Area and perimeter follow
the airway shape (rectangle, round; other profiles are scaled) and `Resistance Total` is
scaled with the Atkinson relation `R ~ k * P * L / A^3`, so shock and density factors Ventsim
folded into it are kept in proportion. `add_airways` copies an existing airway to new end
points; it does not copy fans, fixed flows or pressures, heat or gas sources. A new end
joins the network where it meets an existing airway end (not at a surface connection,
which it would close); an end that meets nothing is a dead end unless listed in
`open_ends`. An airway whose two ends meet the same junctions as an existing one is
refused: Ventsim rejects it as a duplicate. For twin drives set `Pathways` (Ventsim's
Parallel Airways) on the existing airway instead; resistance is scaled by (old/new)^2.
Run `simulate_airflow` on the copy for the new airflow.

## Model checks (Ventsim warnings)

`model_check` lists the Ventsim warnings and errors (manual, Appendix E) that can be found
from the file before opening it: loose ends without a surface or closed-end setting
(no entry / no exit), duplicate airways, zero-length airways, a running fan or fixed flow
whose air has no path once closed airways are taken out (operating fan blocked / fix
stopped), fans in series, stoppings in series, and conflicting fixed flows on one path.
`fix_model` writes a copy with the safe fixes: duplicates merged into one airway with
`Pathways` and the parallel resistance, air-less loose ends marked Closed End, blocked fans
turned off. `edit_airways` reports issues its edits create (closing an airway that is the
only path of a fan, for example).

Checked against Ventsim 6.0.5: the copy with an airway drawn on top of another gives
Ventsim's "2 duplicate airway" and `model_check` reports the same pair; closing the intake
of Tutorial2 gives Ventsim's "operating fan blocked" and `model_check` reports that fan.
On the 20 shipped examples it reports no errors; warnings only in HVACDemoB (loose duct
ends) and RealMine (stoppings in series). Loose ends that carried air in the saved results
are not reported: Ventsim treated those as open. Warnings that need a simulation (stalled
or reversed fans, recirculation, density, heat) come from `fan_duty_check` and
`recirculation_check`.

## Checked in Ventsim

Copies written by this server, re-simulated in Ventsim DESIGN 6.0.5 (airflow) and compared
airway by airway (sum of absolute differences / total airflow):

| Copy | Change | Prediction vs Ventsim |
|---|---|---|
| Tutorial2 (`write_results`) | two resistances x3 and x0.5 | 0.04 % (`what_if`) |
| Tutorial2 (`edit_airways`) | a length and a width changed | 0.20 % (`simulate_airflow`) |
| Tutorial2 (`edit_airways`) | intake closed | both: no airflow anywhere |
| Fire (`write_results`) | fan at 90 % | 9.5 % (`what_if`) |

Ventsim re-simulating Fire unchanged reproduces its saved file within 0.28 %, so the Fire
gap comes from the change: with the fan at 90 % Ventsim also recomputes the fire and rock
heat (inlet temperatures drop by up to 3 K) and with them the natural ventilation, which
`what_if` keeps frozen. This comparison also showed that Ventsim puts fans whose airway
equation uses total pressure on the curve's total-pressure column; `what_if` now does the
same (it scaled the static curve before, 10.7 %).

## Limits

- Fire dynamics, dispersion in time-dependent contaminant spread and rock heat from rock properties stay in Ventsim.
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
