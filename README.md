# Li-ion Battery SOC Estimation: 1-RC ECM + Extended Kalman Filter vs. Coulomb Counting

State-of-charge (SOC) estimation for a real Panasonic NCR18650PF cell,
comparing an Extended Kalman Filter (EKF) built on a 1-RC Equivalent
Circuit Model (ECM) against plain Coulomb counting (CC), across five
conditions: nominal, wrong initial SOC, corrupted sensors, a mismatched
process model, and a wrong assumed capacity. All real laboratory data
throughout (Kollmeyer/Mendeley Panasonic 18650PF dataset, 25°C), no
synthetic cell data.

**Report:** [`report/report.pdf`](report/report.pdf) - start here. It
covers the data, the model, the estimator, the five experiments, and the
main finding: the EKF is not a uniform improvement over Coulomb counting.
It wins when the fault is something CC structurally can't recover from (a
wrong initial SOC), and loses to CC's simplicity when the fault is small
relative to the ECM's own ~72 mV local voltage-fit bias.

## Layout

```
src/            Core, tested functions (data loading, OCV fit, ECM,
                Coulomb counting, EKF, metrics). Each has an arguments
                block, explicit error identifiers, and a matching test
                in tests/.
experiments/    expA_nominal.m ... expE_capacity.m - the five evaluation
                scripts. Each is self-contained: run it to reproduce its
                numbers and figure.
simulink/       ecm_ekf.slx (built by buildModel.m) - the EKF running as
                a Simulink MATLAB Function block, validated against
                src/ekfEstimate.m by runSimulinkValidation.m.
tests/          MATLAB unit tests (testXxx.m) and evidence scripts
                (checkXxx.m) that back specific numbers used in
                docs/DECISIONS.md.
docs/           DECISIONS.md is the project's decision log: what was
                chosen, why, and the evidence behind each number. Read
                this for the reasoning the report only summarizes.
report/         report.md / report.pdf - the final write-up.
figures/        The plots embedded in the report, one per experiment.
scripts/        crosscheck.py - an independent Python re-implementation of
                the same equations, used to verify the MATLAB numbers.
```

## Running it

1. `src/loadDataset.m` onward assume the raw `.mat` files are in
   `data/raw/` (already included).
2. Run any `experiments/expX_*.m` script directly. Each one adds `src/`
   to the path itself, loads the data, runs both estimators, prints
   metrics and an interpretation, and saves a figure to `figures/`.
3. `runSimulinkValidation.m` builds and validates the Simulink model
   (it calls `buildModel.m` itself). Confirmed working, run end-to-end in
   real MATLAB/Simulink on 2026-09-18: exact match,
   `max |SOC_simulink - SOC_ekfEstimate| = 0.000e+00`.
4. All `tests/testXxx.m` files run under MATLAB's `runtests`.
5. `python scripts/crosscheck.py`, run from the repository root, reproduces
   every number in the table below without using any of the MATLAB code.
   Needs numpy and scipy, nothing else.

## Status

Every number below was produced by running the experiment scripts in
MATLAB, and independently reproduced by `scripts/crosscheck.py`, a
separate Python implementation of the same equations that shares no code
with the MATLAB pipeline - two implementations, so a transcription error
in either would show up as a disagreement. Nine of the ten figures agree
exactly to four decimal places. The tenth, Experiment C's EKF value, is
the one case that cannot be compared exactly: its voltage-noise draw
comes from MATLAB's `rng(42)` stream, which numpy does not reproduce
bit-for-bit, and the two implementations land 0.0005 percentage points
apart. Simulink matches `src/ekfEstimate.m` exactly.

| Experiment | Condition | EKF RMSE | CC RMSE | More accurate |
|---|---|---|---|---|
| A | nominal | 7.4317% | 0% (matched, uninformative) | - |
| B | wrong initial SOC | 7.3955% | 20.0% (constant, uncorrectable) | EKF |
| C | current bias + voltage noise | 7.4455% | 0.5157% | CC |
| D | mismatched process model | 8.9161% | 0% (reference point) | - |
| E | capacity mismatch | 7.7436% | 1.7294% | CC |

The figures in `report/report.pdf` and in `figures/` are the plots from
these runs.
