# PIDA GE Evaluation on the Trans-African Network — Workflow

This document records the end-to-end workflow that produced the
General-Equilibrium (GE) evaluation of the PIDA road pipeline in
§7.1 of the paper. It is the trans-African-network counterpart of
the partial-equilibrium analysis in `PIDA_PE_analysis.R`.

> **TL;DR.** We restrict attention to the strategic trans-African
> corridor system (the 47-largest-port-cities fastest-routes graph,
> 330 condensed edges, 212 nodes). Of those 330 edges, 84 (~25 %)
> coincide with a PIDA project. We then solve the Fajgelbaum–Schaal
> welfare-maximising planner at the PIDA-equivalent budget
> $K_i = 10.089$ B USD'15 and overlay the PIDA-flagged subset of
> trans-African links (dotted navy) on the four GE-optimal
> allocations (yellow→red intensity). Frictions are omitted; the
> reason is documented in §4 of the paper (their effect on
> allocation shape is minor).

---

## 1. Why the trans-African network (not the full grid)?

An earlier draft compared PIDA to the GE-optimal allocation on the
full continental grid (12,092 cells, 2,825 edges). That comparison
was misleading: it implicitly assumed the planner could freely
allocate spending across the full grid, while PIDA is restricted by
design to a small set of strategic continental corridors. The
resulting "PIDA captures 2–6 % of optimal budget" headline conflated
the *spatial alignment* question with the much-larger *network
scope* question.

The fair benchmark is to put PIDA and the GE planner on the **same
restricted set of strategic corridors**. The "fastest routes between
the 47 largest port-cities" graph constructed in
`code/9_major_transport_routes.R` is exactly that set: 330 condensed
edges, each representing a real-route segment of the strategic
trans-African network. Within this network, the comparison is
apples-to-apples.

## 2. PIDA → trans-African edge matching

The matching was done in `code/9_major_transport_routes.R`, section
**"Match to PIDA"** (L214–238):

```r
real_paths_simplified <- qs2::qs_read(
  "data/transport_network/trans_african/trans_africa_network_47_largest_fastest_real_edges.qs2"
)
PIDA_bridge   <- fread("data/PIDA/PIDA_edges_bridge.csv")
PIDA_edges    <- join(edges_all_param, PIDA_bridge) |> subset(PIDA == "Yes")
ind           <- st_nearest_feature(st_centroid(PIDA_edges),
                                    st_centroid(real_paths_simplified)) |> unique()
# -> 84 distinct trans-African edges, indices in [1, 330]
```

The 84 indices were pasted verbatim into
`PIDA_GE_analysis_trans_african.R` (variable `PIDA_ind`) so that the
GE script is fully self-contained and does not need to re-run the
geographic matching.

## 3. Julia GE simulations

Run by `code/11_GE_simulation_trans_african/optimal_trans_african_networks_largest_pcities_dual_loop.jl`
with `Ki = 10.089e3` (millions USD'15 → \$10.089 B).
The loop iterates over

| Parameter      | Values used in this evaluation |
| -------------- | ------------------------------ |
| `Ki`           | `10.089e3`                     |
| `gamma`        | `gamma_DRS = 0.946` (standard), `gamma_IRS = β²/γ` (IRS, weak congestion) |
| `rho`          | `0` (utilitarian), `2` (inequality-averse) |
| `sigma`        | `3.8` (Armington baseline; we use this in the paper) |
| `beta`         | `1` (fixed) |
| `cross_good_congestion` | `true` |
| `duality`      | `true` |

The four runs used in the paper are the four $(\gamma,\rho)$
combinations at $\sigma = 3.8$.

### Output filenames

```
results/transport_network/GE_dual/trans_african/
├── nodes_results_22g_10089m_fixed_cgc_sigma3.8_rho0_duality_julia.csv   # standard, ρ=0
├── nodes_results_22g_10089m_fixed_cgc_sigma3.8_rho2_duality_julia.csv   # standard, ρ=2 (IA)
├── nodes_results_22g_10089m_fixed_cgc_irs_na_sigma3.8_rho0_duality_julia.csv  # IRS, ρ=0
├── nodes_results_22g_10089m_fixed_cgc_irs_na_sigma3.8_rho2_duality_julia.csv  # IRS, ρ=2 (IA)
└── edges_results_22g_10089m_fixed_cgc_*.csv                              (4 matching files)
```

The `22g` prefix records the 22 differentiated goods, `10089m`
records the budget in millions, `cgc` records cross-good congestion
on, `irs_na` marks the IRS-without-annealing runs.

## 4. R analysis script

`code/12_PIDA_evaluation/PIDA_GE_analysis_trans_african.R` does the
following for each of the four specs:

1. Loads the GE result CSVs (nodes + edges, 212 × N and 330 rows).
2. Joins the edges to the simplified real-paths geometry
   (`trans_africa_network_47_largest_fastest_real_edges.qs2`) by
   row position — both have 330 rows in the same order.
3. Computes per-edge $\Delta\%$-upgrade
   `perc_ug = (Ijk - Ijk_orig) / (100 - Ijk_orig) * 100`, capped to
   `[0, 100]`.
4. Computes statistics: total budget spent, kilometres worked,
   PIDA-overlap share (by both budget and km), aggregate welfare
   gain ($\sum u_j L_j$ ratio), consumption gain ($\sum C_j$
   ratio), MA gain on the GE-optimised network.
5. Produces two figures per spec:
   - `trans_africa_network_GE_<spec>_perc_ug_with_PIDA.pdf`
     (intensive-margin investment map, PIDA-overlap edges
     overlaid as dotted navy lines)
   - `trans_africa_network_GE_<spec>_upw_gain.pdf`
     (per-node welfare gain, blue→red diverging palette)

The script also writes a summary `qs2` object:

```
results/transport_network/PE/PIDA_GE_trans_african_results.qs2
```

containing `stats_tbl` (one row per spec), `all_stats` (the same
data as a list), the `PIDA_ind` vector, and the network-size
metadata.

## 5. Headline numbers from the run

### σ = 3.8 (Armington baseline; main paper figure)

|             | spec                              | Budget | km    | PIDA-share (\$) | PIDA-share (km) | WG    | CG    | MA    |
| ----------- | --------------------------------- | -----: | ----: | --------------: | --------------: | ----: | ----: | ----: |
| Standard ρ=0 | `..._sigma3.8_rho0_duality_julia` | 10.089 B | 40,889 | **21.2 %**      | 22.5 %          | 0.085 % | 0.10 % | 19.8 % |
| Standard ρ=2 | `..._sigma3.8_rho2_duality_julia` | 10.089 B | 39,459 | **20.4 %**      | 22.4 %          | 0.045 % | 0.06 % | 20.7 % |
| IRS ρ=0      | `..._irs_na_sigma3.8_rho0_*`      | 10.089 B | 40,676 | **21.5 %**      | 22.8 %          | **0.40 %** | 0.48 % | 20.4 % |
| IRS ρ=2      | `..._irs_na_sigma3.8_rho2_*`      | 10.089 B | 39,995 | **21.5 %**      | 23.0 %          | 0.18 % | 0.26 % | **21.3 %** |

Edge-share floor: 84 / 330 ≈ 25.5 % of the network is PIDA-flagged,
so the planner's 20–22 % budget share on PIDA edges is
near-proportional to PIDA's edge count.

### σ = 2.0 (robustness; appendix figure)

|             | spec                              | Budget | km    | PIDA-share (\$) | PIDA-share (km) | WG    | CG    | MA    |
| ----------- | --------------------------------- | -----: | ----: | --------------: | --------------: | ----: | ----: | ----: |
| Standard ρ=0 | `..._sigma2.0_rho0_duality_julia` | 10.089 B | 39,105 | **25.6 %**      | 27.4 %          | 1.11 % | 1.52 % | 20.8 % |
| Standard ρ=2 | `..._sigma2.0_rho2_duality_julia` | 10.089 B | 38,823 | **26.3 %**      | 28.2 %          | 0.69 % | 0.99 % | 21.1 % |
| IRS ρ=0      | `..._irs_na_sigma2.0_rho0_*`      | 10.089 B | 38,560 | **25.8 %**      | 27.7 %          | **2.78 %** | 3.87 % | 20.1 % |
| IRS ρ=2      | `..._irs_na_sigma2.0_rho2_*`      | 10.089 B | 38,005 | **27.7 %**      | 30.2 %          | 2.40 % | 3.45 % | 19.9 % |

Lower σ makes goods less substitutable: the planner spreads spending
more evenly across the strategic network, raising the PIDA-overlap
share to 25.6–27.7 % (about proportional to PIDA's edge count) and
delivering much larger aggregate welfare gains (0.69 %–2.78 % vs
0.04 %–0.40 % under σ = 3.8). MA gains remain ≈20 % in both.

## 6. Where the results appear in the paper

- **Subsection** `\subsection{General Equilibrium Evaluation}`
  (label `sec:GE_PIDA`) of the
  `Evaluating PIDA Road Projects` section.
- **Figure** `\ref{fig:PIDA_GE_main}` — 2×2 panel of the four
  `perc_ug_with_PIDA.pdf` maps.
- Headline stats embedded directly in the panel column headers
  and discussed in the surrounding prose.

## 7. Reproduce

Pre-requisite: the four CSV files in
`results/transport_network/GE_dual/trans_african/` produced by
`optimal_trans_african_networks_largest_pcities_dual_loop.jl`
(already present in the repo).

```bash
cd OptimalAfricanRoads
Rscript code/12_PIDA_evaluation/PIDA_GE_analysis_trans_african.R
```

This will write all 8 figures and the
`PIDA_GE_trans_african_results.qs2` summary. Recompile the paper
with the usual two-pass `pdflatex` to ingest the new figures.

## 8. Decisions / non-trivial choices

1. **Network scope.** The 47-largest-port-cities fastest-routes
   graph was chosen as the GE comparison set, because PIDA is itself
   a continental-corridor pipeline. The grid-level comparison
   inherited from an earlier draft was discarded.
2. **Budget = \$10.089 B.** Matches the PE PIDA cost
   (`sum(pfirst(ug_cost_km, cost_km) * distance / 1000)` over
   `PIDA == "Yes"`, computed in
   `code/12_PIDA_evaluation/PIDA_edges_osbp_export.R`).
3. **Frictions excluded.** Justified by §4 of the paper, which
   shows frictions move levels but not the spatial allocation
   shape; including them would not change the PIDA-vs-optimal
   comparison.
4. **OSBP excluded.** Their effect is already cleanly identified in
   PE (the +3.6 pp MA gain from halving border friction at the 27
   OSBP candidate edges); the GE counterfactual focuses on the
   road-investment allocation only.
5. **rho=2 utility re-mapping.** When `rho == 2` the model
   solver returns `uj` as the negative reciprocal of true utility;
   the script applies `uj <- (uj * (-1))^(-1)` before computing
   ratios (mirrors `analyze_trans_african_results.R` L95).
6. **MA gain semantics.** Computed on the GE-optimised network as
   the ratio of total-MA at the new network's shortest-path
   travel-times divided by total-MA at the original. This is the
   same definition used throughout the PE sections and in
   `analyze_trans_african_results.R`.
