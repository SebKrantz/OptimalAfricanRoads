# ATRC-2027 — Form Answers

**Paper:** Optimal Investments in Africa's Road Network
**Author:** Sebastian Krantz (Kiel Institute for the World Economy and World Bank)
**Selected theme:** 8. Data, Analytical Tools and Emerging Technologies (with relevance to 6. Policy, Planning and Governance)

## Section 2 — Research background and aims

Africa's road network is sparse, fragmented, and burdened by some of the world's highest cross-border trading costs. Yet over US$160 billion of road and border-post investments are on the table, most prominently through the African Union's Programme for Infrastructure Development in Africa (PIDA). Existing studies evaluate selected corridors or planned packages, but none characterise what an economically optimal trans-African road investment programme would actually look like, nor benchmark PIDA against such an optimum. This paper fills that gap. It aims to (i) take a comprehensive empirical stocktake of Africa's present road network; (ii) characterise globally optimal road investments from both market-access and welfare perspectives; (iii) quantify how cross-border frictions reshape optimal spatial allocations; and (iv) use these benchmarks to evaluate the PIDA pipeline of upgrades and one-stop border posts (OSBPs).

## Section 3 — Research method

The empirical backbone is 144 million OSRM-derived trans-continental shortest-path routes computed on a topographically detailed road graph, used to measure local network efficiency and continental market access (MA). A condensed graph connects 447 cities (population >100,000) and 52 major ports through 1,379 nodes and 2,344 edges; an algorithm proposes 481 candidate new links. Cross-border friction data come from Doing Business and AUDA-NEPAD OSBP candidate locations; PIDA projects are encoded as edges. Two complementary analytical frameworks are applied. The partial-equilibrium framework simulates MA-maximising upgrades, ranks links by dollar-per-minute-of-MA-per-dollar ($/min/$), and conducts macroeconomic cost-benefit analysis. The general-equilibrium framework adapts Fajgelbaum and Schaal's (2020) optimal-transport-network model with an Armington elasticity of 3.8 and 22 differentiated goods, run under both standard and increasing-returns-to-infrastructure (IRS) planners and under utilitarian and inequality-averse social welfare functions. The same machinery evaluates PIDA against cost-matched optimal counterfactuals.

## Section 4 — Preliminary findings

Cross-border frictions are the dominant constraint on continental connectivity, decisively reshaping where road investments are most valuable. Upgrading existing links yields large MA returns, while new construction is generally lower-yield. PIDA's 187 upgrades and 7 new links cost ~$10.1B and yield an 8.5% frictionless MA gain (14.1 $/min/$); halving border friction at all 27 OSBP-adjacent edges adds 3.6 percentage points. A cost-matched optimal package delivers 18.4%—more than double PIDA's. In general equilibrium, a welfare-maximising planner with the same budget allocates only 2.1%–5.7% to PIDA-flagged links. Macroeconomic break-even requires just 0.06%–0.09% additional growth.

## Section 5 — Stage of research completion

The paper is in advanced revision. A near-complete ~70-page draft exists; PE and GE analyses are finalised, and the PIDA evaluation section was added in May 2026. Reproducibility code, data, and a Julia implementation of the Fajgelbaum–Schaal toolbox are public on GitHub. Remaining work is journal-revision polishing and exposition. The full-paper deadline of 30 November 2026 will be met comfortably.

---
*Word counts (excluding section headings; limits in parentheses):*
- Section 2: 130 words  (≤ 200)
- Section 3: 138 words  (≤ 300)
- Section 4: 92 words   (≤ 100)
- Section 5: 60 words
- **Total: 420 words**
