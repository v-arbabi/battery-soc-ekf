# Decision Log

The decisions that shape this project and that need to hold up in the
report or under questioning. Fifteen entries, in the order the project
encountered them.

Format per entry: **Decision, Reason, Evidence, Affected files.**

---

## 1. Data, conventions, and sign

**Decision:** Use the Panasonic 18650PF nominal-25°C C/20 file (source
rows 1–2452 of 2453; row 2453 is a 10°C reading and is excluded) and the
US06 drive-cycle file, both loaded only through `src/loadDataset.m`.
Normalise current sign exactly once inside that function: `I = -meas.Current`,
so that **positive current means discharge**. Remove only the second row
of exact zero-`dt` duplicate pairs; reject, rather than silently resolve,
any non-identical zero-`dt` pair.

**Reason:** Normalising in one tested place means every downstream
function can assume the sign is right instead of re-deriving it, so a sign
error cannot reappear somewhere else. Discharge-positive also matches the
convention used in the reference textbook equations, which makes the
implementation easier to check against the literature.

On the zero-`dt` pairs: when two rows carry the same timestamp but
different values, there is no correct way to choose between them, and any
choice is a guess. Code that guesses and continues hides its error in the
final result. Code that stops forces the problem into the open.

**Evidence:** `tests/testLoadDataset.m` — both datasets load with strictly
increasing time and the documented row counts.

**Affected files:** `src/loadDataset.m`, `tests/testLoadDataset.m`.

---

## 2. Capacity and SOC anchors

**Decision:** Use the **measured** discharge capacity `Q_Ah = 2.997393193 Ah`
(full-branch right-endpoint trapezoidal integration on the C/20 file,
confirmed within 0.01% by the dataset's own `Ah` counter) as the capacity
for all state equations and reference-SOC construction. The documented
nominal `2.9 Ah` is kept separately, as an "assumed vs measured" mismatch
condition for Coulomb-counting robustness evidence, which feeds
Experiment E. `SOC_ref(1) = 1` on the US06 record is a
**protocol-supported working assumption**, not a directly observed anchor:
the record starts mid-protocol with no observed rest, and the report says
so explicitly instead of treating it as ground truth.

**Reason:** Measured capacity on the same cell is stronger evidence than a
datasheet nominal value. Labelling the SOC(1)=1 assumption as an
assumption protects the project from an unsupported certainty claim — and
that caveat applies to every RMSE and bias number in the report, since all
of them are measured against `SOC_ref`.

**Evidence:** `tests/checkUS06ReferenceSOC.m` (2.586 Ah net discharge,
+0.0056% vs dataset `Ah`), `src/referenceSOC.m`.

**Affected files:** `src/referenceSOC.m`, `src/coulombCount.m`,
`src/ecmParams.m`, `src/ecmParamsHPPC.m`.

---

## 3. OCV(SOC) model

**Decision:** Build the OCV-SOC target as 101 points. Interior points
(SOC 0.01–0.99) are the arithmetic mean of independently normalised
discharge and charge C/20 branches, labelled *low-rate charge–discharge
voltage separation*, not equilibrium hysteresis — hysteresis itself stays
unmodeled. The SOC=0 and SOC=1 anchors use only rested, not loaded,
voltage readings. Fit with a degree-9 polynomial, and take `dOCV/dSOC`
analytically via `polyder` instead of numerically differentiating the raw
target.

**Reason:** A polynomial fit gives a smooth, analytically differentiable
curve. Degree 9 is the lowest degree that came out both monotonic and
inside the fit-error bar; degree 8 was non-monotonic, which is physically
meaningless here — adding charge should not lower the equilibrium voltage —
and can push the filter in the wrong direction. The EKF's measurement
Jacobian needs a smooth derivative, and a numerically differentiated raw
target introduces steps that can destabilise the filter.

The endpoints are treated differently on purpose: at SOC 0 and 1 the cell
genuinely comes to rest and a rested voltage is available, so using a
loaded reading there would shift the anchors of the whole curve.

**Evidence:** `tests/testBuildOCV.m` — RMSE 6.7 mV, max target error
22.3 mV, minimum derivative 0.354 V/SOC over [0,1], all positive, so
monotonic.

**Affected files:** `src/buildOCV.m`, `src/buildOCVTarget.m`.

---

## 4. Two ECM parameter sets, and which is nominal

**Decision:** Two 1-RC parameter sets are kept, on purpose:

- `src/ecmParams.m`: literature-derived (same-cell, 25°C, SOC=0.5,
  FFRLS/UDDS table): `R0=4.500052 mΩ`, `R1=28.981830 mΩ`, `C1=218.873 F`.
- `src/ecmParamsHPPC.m`: identified from the local 5-pulse HPPC file
  through a bounded, capped, three-gate procedure — qualify windows, then
  one global least-squares pass over 64 clean pulses, then one validation
  pass: `R0=31.2245 mΩ`, `R1=11.3288 mΩ`, `tau=5.2396 s`, `C1=462.505 F`.

Local US06 forward-voltage fit (sample-wise, `V_model - V_measured`):

| Parameter set | bias | RMSE | max abs error |
|---|---|---|---|
| Literature | +89.2 mV | 114.3 mV | 688.9 mV |
| HPPC-identified | +71.7 mV | 83.8 mV | 575.3 mV |

**The HPPC-identified set is the nominal parameter set** for
`ekfEstimate.m` and Experiments A, B, C, E. **The literature set stays
unchanged and is reused as the mismatched model for Experiment D**, which
uses a real, published alternative instead of an artificial parameter
perturbation.

**Reason:** Neither set removes the residual bias; both are real and
sourced. Using the better-fitting one as nominal and the worse one as the
"wrong model" condition turns unavoidable model imperfection into
legitimate experimental content instead of a blocker.

**Evidence:** `tests/checkUS06HPPC1RCVoltage.m`,
`tests/checkUS06ECMVoltage.m`.

**Constraint carried forward:** no further HPPC re-identification, window
search, or optimiser re-run. If a third parameter set is ever needed, it
has to come from new same-cell evidence, not a manual tweak. Running an
optimiser repeatedly until it produces a preferred number is no longer
identification; fixing this constraint before seeing the result is what
keeps the distinction defensible.

**Affected files:** `src/ecmParams.m`, `src/ecmParamsHPPC.m`,
`src/estimateHPPC1RC.m`.

---

## 5. Right-endpoint timing convention

**Decision:** For every recurrence in the project (`referenceSOC`,
`coulombCount`, `ecmSimulate`, `ekfEstimate`), sample `k≥2` uses the
interval `[t(k-1), t(k)]` and current `I(k)`, so state, RC branch, and
model voltage at index `k` are all computed on the same interval with the
same current sample. `V1(1) = 0` by convention; the first
model-voltage / first-innovation sample is `NaN` / undefined and is
excluded from all fits and metrics.

**Reason:** Mixing left- and right-endpoint conventions between functions
introduces a one-sample lag that looks like parameter error instead of a
timing error — which means hours spent chasing the wrong quantity. Making
every function agree removes that failure mode by construction.

**Evidence:** `tests/testEcmSimulate.m` and `tests/testEkfEstimate.m` both
assert the right-endpoint recurrence directly.

**Affected files:** `src/referenceSOC.m`, `src/coulombCount.m`,
`src/ecmSimulate.m`, `src/ekfEstimate.m`.

---

## 6. Three SOC trajectories stay structurally separate

**Decision:** `SOC_ref`, `SOC_CC`, and `SOC_EKF` stay three separate
variables, built by three functions with disjoint input contracts.
`referenceSOC` and `coulombCount` never receive each other's output,
voltage, OCV, or ECM quantities, and `ekfEstimate` never receives
`SOC_ref`. Each robustness experiment (B–E) builds `SOC_ref` from clean
inputs and `SOC_CC`/`SOC_EKF` from the corrupted inputs, in that order, so
neither estimator can see the answer.

**Reason:** This is the methodological integrity of the project. If any
function could read the reference it is being scored against, every
comparison in the report becomes circular.

The guarantee is structural, not procedural: `coulombCount` has no
parameter through which voltage or OCV could reach it, and `ekfEstimate`
has no parameter through which `SOC_ref` could. There is nothing to
remember to avoid.

**Evidence:** every `check*.m` / `exp*.m` script's function-call ordering,
enforced by each function's own input list.

**Affected files:** `src/referenceSOC.m`, `src/coulombCount.m`,
`src/ekfEstimate.m`, all `experiments/exp*.m`.

---

## 7. Starting the EKF on an imperfect ECM

**Decision:** Do not require the ECM to pass a voltage-adequacy threshold
before the EKF is written. Proceed with `ecmParamsHPPC()` as the nominal
model (see #4).

**Reason:** The premise of that requirement is backwards. The purpose of a
Kalman filter is to produce a usable estimate *despite* an imperfect
process model; a biased model is a normal input to it, not a precondition
that must be fixed first. No real BMS has a perfect model, and a project
that waits for one never writes a filter.

Practically, both parameter sets in hand are real, sourced, and within a
similar range of residual bias. A third, better set would have required
new same-cell measurements, not another pass over the same data.

That decision is also what produced the project's main finding: because
the filter was allowed to run on a model with a known bias, the
relationship between model fidelity and achievable SOC accuracy became
measurable instead of hypothetical.

**Evidence:** table in #4.

**Affected files:** `src/ekfEstimate.m`, `experiments/expA_nominal.m`.

---

## 8. Experiment A nominal result: read this before the report's abstract

**Finding:** With the nominal (HPPC) parameters, correct initial SOC, and
no added sensor corruption, the EKF's SOC estimate carries a **systematic
bias of roughly −7%** (RMSE ≈ 7.4%, max abs error ≈ 12.6%) relative to
`SOC_ref`. The ±3σ uncertainty band covers under 1% of samples, well
short of the expected ~99.7%.

This is not a filter defect. It was checked against the forward-voltage
residual from #4: the HPPC ECM's own systematic voltage bias (+71.7 mV),
divided by the mean `dOCV/dSOC` over the visited SOC range (≈0.90 V/SOC),
predicts almost exactly this SOC bias (≈8%). The EKF has no way to know
the voltage residual comes from an imperfect process model instead of a
wrong SOC, so it does the only thing available to it: it moves SOC until
the voltage prediction matches. The near-zero ±3σ coverage is the direct
symptom of that. `Q` was set assuming random measurement and process
noise, not a large deterministic model bias, so the filter ends up
overconfident.

`Q` was not inflated to widen the band. A wider band bought by inflating
`Q` after the bias was already visible would be tuning to the answer, and
would trade a traceable miscalibration for one that merely looks better.

That the error is one-sided rather than scattered is itself diagnostic:
bias (−6.8%) and RMSE (7.4%) are close, so almost all of the error comes
from a single systematic offset with a locatable cause.

**Report implication:** Do not claim "EKF RMSE 1–2% in the nominal case."
The real nominal-case finding is that *EKF SOC accuracy is bounded by
process-model fidelity: with the best available same-cell parameter set,
a ~72 mV residual ECM voltage bias — plausibly unmodeled OCV hysteresis
and SOC/temperature-dependent R0/R1, both out of scope, see #9 — converts
to a comparable persistent SOC bias.* It does not change the comparative
story for Experiments B–E, where `SOC_CC`'s input is additionally
corrupted while the EKF fights the same fixed model bias plus that
corruption. That comparison is the content of the project.

**Evidence:** `experiments/expA_nominal.m` output; cross-check against
`tests/checkUS06HPPC1RCVoltage.m` and `tests/checkUS06ECMVoltage.m`
residual statistics.

**Affected files:** `experiments/expA_nominal.m`, `report/`.

---

## 9. Scope exclusions

Not attempted, and stated in the report as scope:

- **Temperature dependence.** Everything is 25°C only. ECM parameters vary
  strongly with temperature; covering it means one parameter set per
  temperature plus an interpolation layer.
- **OCV hysteresis as a model term.** Branch-averaged; the ~45 mV mean
  interior charge/discharge separation is measured and reported, not
  modeled inside the ECM. It is one of the two main suspects behind the
  residual voltage bias, so measuring and declaring it is more useful than
  absorbing it into a fitted correction term.
- **Further HPPC-based parameter identification** beyond the capped
  three-gate procedure in #4.
- **SOH estimation.** Experiment E shows *why* capacity error matters for
  SOC accuracy; the project does not build an SOH estimator. In an
  empirical ECM, capacity fade enters only as an imposed parameter, never
  as a state the model can resolve.
- **Real-time and embedded implementation.** Fixed-point arithmetic and
  execution-time budgets are not addressed.

**Affected files:** `report/` (limitations section), `README.md`.

---

## 10. Experiment B result: wrong initial SOC

`experiments/expB_initSOC.m`. Both `SOC_CC` and the EKF are given
`SOC0 = 0.8` against the true `SOC0 = 1.0` (20 percentage points off);
current and voltage are otherwise clean. `tuning.P0_SOC = (0.20)^2`,
larger than Experiment A's `(0.01)^2`, because here the initial guess is
known to be untrustworthy, and a small `P0_SOC` would make the filter
trust a wrong prior too much.

That change is set from what is known about the experiment before it runs,
not from its output — which is what separates it from tuning to the
answer.

**Result:**

- CC: constant −0.2 error for the entire record; initial error equals
  final error to numerical precision. Coulomb counting has no mechanism to
  correct an initial-condition error.
- EKF: starts at the same −0.2 error and is pulled back by the voltage
  update over time. Final error ≈ −0.084, mean |error| over the back half
  of the record ≈ 0.093. It converges toward Experiment A's steady-state
  bias level (≈ −0.068), not to zero, because the same ~72 mV HPPC-ECM
  voltage bias is still present.
- This is the clean illustration of why an estimator is worth having over
  Coulomb counting: it recovers from a bad initial guess, and CC
  structurally cannot.

**Affected files:** `experiments/expB_initSOC.m`.

---

## 11. Experiment C result: corrupted sensors

`experiments/expC_sensor.m`. A constant `+0.02 A` current bias, plus
zero-mean Gaussian voltage noise (3 mV std, matching the EKF's own `R_V`),
are fed to both `SOC_CC` and the EKF; `SOC_ref` stays clean.

**Result:**

- CC: RMSE ≈ 0.52%, max error ≈ 0.89%. Small, because a `+0.02 A` bias is
  a small fraction of the drive-cycle current.
- EKF: RMSE ≈ 7.4%, almost identical to Experiment A's nominal case.
- **CC is more accurate than the EKF in this condition.** The EKF's error
  is dominated by the pre-existing ~72 mV HPPC-ECM voltage-fit bias from
  Experiment A, not by the sensor fault added here, so fusing a
  noisy/biased current signal with voltage does not help: the voltage
  channel already carries a larger, unrelated error source.
- Reported as measured: a case where sensor fusion does not beat the
  simple baseline, and the reason why.

The script's interpretation text and plot title were written before the
run and assumed the opposite outcome. They were corrected to match the
measured result.

**Affected files:** `experiments/expC_sensor.m`.

---

## 12. Experiment D result: mismatched process model

`experiments/expD_modelMismatch.m`. The EKF is given the
literature-derived constant parameter set (`src/ecmParams.m`, SOC = 0.5
constants) instead of the HPPC-identified nominal set; current, voltage,
and initial SOC are otherwise identical to Experiment A.

**Result:**

- EKF: RMSE ≈ 8.9%, bias ≈ −8.0%, max error ≈ 16.2%, worse than
  Experiment A's 7.4% RMSE and −6.8% bias. Consistent with the literature
  set's larger local US06 voltage-fit error (+0.089 V bias, 0.114 V RMSE
  against the HPPC set's +0.072 V and 0.084 V; see
  `tests/checkUS06ECMVoltage.m` and `tests/checkUS06HPPC1RCVoltage.m`).
- CC is unaffected — same clean current, same capacity as Experiment A —
  and is shown only as a fixed reference point, not as a comparison target.
- This is the cleanest of the five results: a larger model-voltage error
  measurably produces a larger EKF SOC error, isolated from any other
  confound.

**Affected files:** `experiments/expD_modelMismatch.m`.

---

## 13. Experiment E result: capacity mismatch

`experiments/expE_capacity.m`. Both `SOC_CC` and the EKF are given the
documented nominal capacity (`2.9 Ah`) instead of the measured capacity
(`2.997393193 Ah`) used to build `SOC_ref`; current, voltage, and initial
SOC are otherwise clean.

**Result:**

- CC: RMSE ≈ 1.73%, bias ≈ −1.49%, the expected throughput-scaled drift
  from a ~3.2% capacity error.
- EKF: RMSE ≈ 7.74%, bias ≈ −7.11%, only slightly worse than Experiment
  A's 7.43% RMSE.
- **CC is again more accurate than the EKF**, for the same reason as
  Experiment C: the EKF's error is still dominated by the ~72 mV HPPC-ECM
  voltage bias present since Experiment A, and the capacity mismatch adds
  only a small increment on top of it.
- The EKF does not correct the capacity error just because it also uses
  voltage. Its own process model still uses the same wrong `Q_As` in the
  SOC prediction step; the voltage update applies a correction to the
  state, but the prediction keeps advancing on the wrong capacity.
  Correcting it would mean estimating capacity as a state, which is the
  SOH problem this project explicitly does not enter (#9).

**Affected files:** `experiments/expE_capacity.m`.

---

## Cross-experiment summary

| Experiment | Condition | EKF RMSE | CC RMSE | More accurate |
|---|---|---|---|---|
| A | nominal | 7.43% | 0% (matched integration, uninformative) | — |
| B | wrong initial SOC | 7.40% | 20.0% (constant, uncorrectable) | EKF |
| C | current bias + voltage noise | 7.44% | 0.52% | CC |
| D | mismatched process model | 8.92% | 0% (unaffected, reference point) | — |
| E | capacity mismatch | 7.74% | 1.73% | CC |

The defensible thesis is **not** "the EKF always wins." It is this: the
EKF recovers from a bad initial condition where CC structurally cannot
(B), but in this project's setup the EKF's accuracy is bounded below by
the ECM's own ~72 mV local voltage-fit bias. When the injected corruption
(C, E) is small relative to that bias, the much simpler CC baseline ends
up more accurate.

---

## 14. Simulink model scope

**Decision:** `simulink/ecm_ekf.slx` (built by `simulink/buildModel.m`) is
one MATLAB Function block running the exact `src/ekfEstimate.m`
per-sample step logic, driven by From Workspace inputs built from the same
US06 record and OCV model as `experiments/expA_nominal.m`.
`simulink/runSimulinkValidation.m` runs it and checks its logged SOC
against `src/ekfEstimate.m`'s output to `1e-9`.

**Reason:** A full block-diagram rebuild — Discrete-Time Integrator,
1-D Lookup Table, Discrete Transfer Fcn — would re-express mathematics
already implemented and tested in MATLAB, in a second notation, with no
new result. The requirement that matters is an estimator running inside
Simulink and numerically validated against the MATLAB reference, at a size
the project can still explain end to end.

**Getting it working took four rounds of correction**, each driven by a
real error message:

1. A fixed-step discrete solver steps at a uniform dt that does not line
   up with the real, nonuniform US06 sample times, so it lands on times
   the input data has no value for. All three From Workspace blocks threw
   *"Unable to extrapolate output values after final workspace data value
   because interpolation is not enabled."*
2. Switching to a variable-step discrete solver ran without error but
   logged only 51 of 48060 samples. A variable-step discrete solver with
   no block explicitly registering a sample hit falls back to an internal
   `MaxStep` heuristic, `(StopTime-StartTime)/50`, which ignores the
   irregular times in the From Workspace data entirely. It was never
   stepping at the real US06 sample times, just about 50 evenly spaced
   guesses.
3. Redesigned to take real time out of Simulink's clock entirely: the
   From Workspace blocks now use the plain sample index (0, 1, …, N−1, one
   integer tick per row) as their time base, with an explicit
   `SampleTime=1` and a fixed-step discrete solver, `FixedStep=1`. That
   guarantees exactly one simulation step per data row. The real elapsed
   time between samples, needed for `alpha=exp(-dt/tau)`, was already
   carried as an explicit *data* value in `dt_ws`'s second column instead
   of derived from Simulink's clock, so this changes only how many steps
   Simulink takes, not the EKF mathematics.
4. That still threw the same extrapolation error, even though the data
   covers exactly `[0, N-1]` and `StopTime = N-1`: Simulink queries at
   least one time point outside that closed interval for termination
   bookkeeping, whatever the solver type. Fixed by enabling
   `Interpolate data` on all three blocks. Since the supplied data already
   has an exact entry at every integer tick, this changes none of the
   values used at 0..N−1; it stops Simulink erroring on what it queries
   just outside that range.

One further environment detail: this MATLAB version returns To Workspace
outputs as fields of the `Simulink.SimulationOutput` object
(`out.SOC_simulink`), not as base-workspace variables.
`runSimulinkValidation.m` reads `simOut.SOC_simulink` accordingly, with a
base-workspace fallback for older releases.

**Confirmed result:**

```
Simulink vs. MATLAB (src/ekfEstimate.m) EKF SOC comparison, Experiment A conditions
  max |SOC_simulink - SOC_ekfEstimate| = 0.000e+00
  PASS: Simulink model matches the MATLAB reference to 1e-09.
```

Bit-for-bit agreement, not merely within tolerance.

**Affected files:** `simulink/buildModel.m`, `simulink/runSimulinkValidation.m`.

---

## 15. All five experiments confirmed in MATLAB

`experiments/expA_nominal.m` through `expE_capacity.m` were run in MATLAB
and the printed console output recorded for each. Every metric matched the
independent Python implementation of the same equations (#8, #10–13) to
rounding:

| Experiment | Metric | Python implementation | MATLAB |
|---|---|---|---|
| A | EKF RMSE / CC RMSE | ≈7.43% / 0% | 7.4317% / 0.0000% |
| B | EKF RMSE / CC RMSE | ≈7.40% / 20.0% | 7.3955% / 20.0000% |
| C | EKF RMSE / CC RMSE | ≈7.44% / 0.52% | 7.4455% / 0.5157% |
| D | EKF RMSE / CC RMSE | ≈8.92% / 0% | 8.9161% / 0.0000% |
| E | EKF RMSE / CC RMSE | ≈7.74% / 1.73% | 7.7436% / 1.7294% |

No discrepancy exceeded rounding.

This matters beyond confirming the figures. The same equations were
implemented twice, independently, in two languages; a transcription error
in either would show up as a disagreement between them. Combined with the
bit-for-bit Simulink agreement above, every technical component of the
project — data handling, ECM, EKF, all five experiments, and the Simulink
implementation — has been checked against at least one independent
implementation.
