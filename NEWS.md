# caabgp 0.2.0

- Align sequential allocation with the manuscript algorithm by requiring the
  entire prespecified cohort to satisfy both sample-size and budget limits.
- Add stratum-specific acquisition thresholds and consecutive-update stopping
  counters, with default patience `J + 1`.
- Update acquisition values and stopping counters for every active stratum
  after each completed cohort.
- Retain terminal recommendations for all strata using the complete trial
  history, including strata closed to further exploration.
- Align the two-agent empirical-Bayes optimizer starts and default iteration
  limit with the manuscript simulation implementation.
- Require the initial sample to include at least two patients per stratum.
- Add regression tests for full-cohort feasibility, stratum closure, and final
  all-stratum recommendations.
