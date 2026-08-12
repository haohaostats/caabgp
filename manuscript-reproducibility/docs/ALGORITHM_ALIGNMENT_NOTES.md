# Algorithm-aligned simulation branch (2026-08-06)

This directory is an independent copy of the empirical-Bayes R pipeline used
for the current manuscript.  The May 2026 pipeline and its results remain
unchanged in their original directory.

## Manuscript-to-code alignment

The sequential implementation follows the current `manuscript/manuscript.tex`
algorithm in three places that were not all represented in the earlier
simulation pipeline:

1. A candidate action is feasible only when the *entire* prespecified cohort
   can be enrolled without exceeding either `N_max` or `B_max`.  The cohort is
   never reduced at the last allocation.
2. Acquisition values and stopping counters are updated for every active
   stratum after every newly observed cohort.  A stratum closes after the
   threshold condition holds at the prespecified number of consecutive interim
   analyses.  The primary simulations retain `delta_k = 0`, so threshold-based
   early stopping is disabled there, as stated in Supplementary Table S1.
3. The terminal fit uses the complete trial history and predictions are
   evaluated for all strata, including strata closed earlier.  Final
   recommendations therefore continue to benefit from the shared GP component.

The pooled S-GP comparator also obeys the strict full-cohort budget and sample
size constraints.  Personalized comparators use the same feasible
stratum--dose action construction as the proposed design, with their respective
surrogate and acquisition definitions.

## Validation performed before the full run

- All R sources loaded without parse errors.
- `validate_algorithm_alignment.R` passed its fixed-cohort, active-stratum,
  stopping-counter, and final-all-stratum checks.
- A two-replicate run across all four scenarios and all five methods completed;
  every recorded adaptive allocation had `cohort_size == planned_cohort_size`
  and no sample-size or budget violation occurred.
- The real-data-calibrated application completed under the aligned code.  Its
  application CSV and LaTeX summaries were byte-identical to the earlier run,
  because all application allocations already used complete cohorts of two.

## Full simulation command

```powershell
& 'D:\Soft\R-4.5.2\bin\x64\Rscript.exe' run_full_simulation_study_R.R 1000 12 TRUE
```

The portable launcher writes a timestamped run log under `logs/`.

For portable Windows use, `RUN_PORTABLE.ps1` and the two `.bat` launchers add
automatic R detection, worker selection, elapsed-time/ETA progress, log capture,
result validation, and worker-aware chunk checkpoints for safe restart.
