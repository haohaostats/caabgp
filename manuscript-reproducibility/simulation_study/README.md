# Simulation study

The four scenario directories follow the manuscript order:

1. Scenario 1: shared main surface with moderate subgroup deviations (Figure 2, Table 3).
2. Scenario 2: rare, noisy, clinically prioritized stratum (Figure 3, Table 4).
3. Scenario 3: mixed heterogeneity, rare-stratum prioritization, and high manufacturing cost (Figure 4, Table 5).
4. Scenario 4: high outcome noise and value of replication (Figure 5, Table 6).

The `sensitivity_analysis/` directory contains the Scenario 3 sensitivity results (Table 7). Each directory includes the verified 1,000-replicate compact numerical summaries and the corresponding publication output.

Scenario definitions and response surfaces are fixed in `R/00_settings.R` and `R/01_surfaces.R`. Run `scripts/01_run_full_simulation.R` to repeat the Monte Carlo study, or `scripts/02_rebuild_main_text_outputs.R` to recreate the displays from the included verified summaries.

