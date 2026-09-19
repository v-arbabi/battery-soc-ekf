# First-Order Thevenin ECM Parameter Evidence Proposal

## Status and scope

This document proposes one provisional, source-backed constant parameter set
for the first-order Thevenin model of a Panasonic NCR18650PF cell near 25 degC.
It is supporting evidence for the parameter choice, not approval to
implement or use the parameters. No parameters were fitted from the local C/20
or US06 files.

## Primary source

The primary source is the open-access Energy Reports article
[Research on pulse charging current of lithium-ion batteries for electric
vehicles in low-temperature environment](https://doi.org/10.1016/j.egyr.2023.04.226),
Energy Reports 9, Supplement 7 (2023), pages 1447-1457.

The article identifies a first-order Thevenin model for the Panasonic
NCR18650PF using voltage and current from the Kollmeyer Panasonic 18650PF UDDS
dataset. It reports that the 25 degC UDDS data have a 0.1 s sampling interval
and uses forgetting-factor recursive least squares (FFRLS) to identify the
model parameters. Its Table 2 reports parameters at SOC values from 0.1 to
0.9. The cell identity and dataset lineage match this project's selected cell
family and source dataset.

The local dataset description independently identifies the tested cell as a
brand-new 2.9 Ah Panasonic 18650PF cell and documents the 25 degC test series:
[`data/raw/Readme file - desc of tests performed.txt`](../data/raw/Readme%20file%20-%20desc%20of%20tests%20performed.txt).

## Source model and parameter mapping

The source uses the following continuous-time first-order Thevenin equations:

```text
dU_D/dt = i_L/C_D - U_D/(C_D*R_D)
U_t     = U_OC - U_D - i_L*R_i
```

The project's intended model is:

```text
dV1/dt = I/C1 - V1/(R1*C1)
V_t    = OCV(SOC) - V1 - I*R0
```

The topology and equations support this direct symbol mapping:

| Project parameter | Source symbol | Definition | Unit | Evidence type |
|---|---|---|---|---|
| `R0` | `R_i` | Instantaneous/ohmic series resistance | ohm | Directly reported; symbol renamed |
| `R1` | `R_D` | Polarisation resistance of the single RC branch | ohm | Directly reported; symbol renamed |
| `C1` | `C_D` | Polarisation capacitance of the single RC branch | farad | Directly reported; symbol renamed |
| `tau` | `R_D*C_D` | RC polarisation time constant | second | Source-defined relation; numerically calculated here |

No unit scaling is needed because the source reports resistance in ohms and
capacitance in farads. For readability only, resistance is also shown below in
milliohms using `1 ohm = 1000 milliohms`.

## Proposed constant parameter set

The proposed fixed set uses the source's 25 degC, `SOC = 0.5` row:

| Parameter | Proposed value | Equivalent display | Provenance |
|---|---:|---:|---|
| `R0` | `0.004500052 ohm` | `4.500052 milliohm` | Direct source value `R_i` at 25 degC and SOC 0.5 |
| `R1` | `0.028981830 ohm` | `28.981830 milliohm` | Direct source value `R_D` at 25 degC and SOC 0.5 |
| `C1` | `218.8730953 F` | - | Direct source value `C_D` at 25 degC and SOC 0.5 |
| `tau` | `6.343342840 s` | - | Calculated from the three source values |

The required calculation is:

```text
tau = R1 * C1
    = 0.028981830 ohm * 218.8730953 F
    = 6.343342839558 s
```

The rounded proposal records `tau = 6.343342840 s` while retaining the direct
source precision for `R1` and `C1`.

## Why the SOC 0.5 row is proposed

The source reports SOC-dependent parameters rather than one universal constant
set. Selecting its `SOC = 0.5` row is an explicit engineering simplification
for this project's constant-parameter 1-RC implementation:

- it is the midpoint of the physical SOC interval;
- it is a directly tabulated same-cell, 25 degC operating point;
- the source states that polarisation behaviour is most favourable around
  50% SOC and worsens toward low and high SOC;
- it avoids inventing an average or fitting a new SOC-dependent parameter
  function from published table values.

The selection of the 50% row is an inference and proposal. The numerical
`R_i`, `R_D`, and `C_D` entries themselves are directly reported by the source.

## Identification conditions and compatibility

| Item | Source condition | Project interpretation |
|---|---|---|
| Cell | Panasonic NCR18650PF, nominal 2.9 Ah | Same cell type and dataset lineage |
| Temperature | 25 degC table | Matches the project's nominal 25 degC scope |
| SOC | Proposed row is SOC 0.5; source table covers 0.1-0.9 | One fixed midpoint set, not an SOC-varying model |
| Excitation data | UDDS voltage and current | Dynamic drive-cycle identification, not local C/20 fitting |
| Identification method | FFRLS | Values are literature-identified, not re-identified here |
| Model order | One RC polarisation branch plus series resistance | Matches the planned first-order Thevenin topology |
| Source sampling | 0.1 s for 25 degC UDDS | Compatible with the dynamic dataset class; no resampling claim |

## Limitations and review cautions

1. The source shows substantial SOC dependence. A fixed 50% parameter set is
   a deliberate simplification and must not be described as valid uniformly
   over all SOC values.
2. The proposed `R0` is lower than the broad 20-40 milliohm order-of-magnitude
   the project's scope. The source is same-cell and
   condition-specific, but this discrepancy should be checked during voltage
   validation rather than hidden or manually corrected.
3. The local repository currently lacks the selected cell's HPPC/pulse file,
   so these values cannot yet be independently reproduced from local pulse
   data.
4. Cell age, exact test sequence, sign convention used internally by the
   FFRLS implementation, and estimator sensitivity can affect identified
   values. The source equations establish compatible voltage-drop signs, but
   do not remove ordinary cross-cell and identification uncertainty.
5. This proposal does not approve SOC-dependent interpolation, temperature
   dependence, parameter optimisation, or a second RC branch.

## Recommendation for the gate

Subject to independent verification against the source, use the
25 degC, SOC 0.5 row above as the initial constant 1-RC parameter proposal.
Treat it as literature-derived and provisional until the later model-voltage
comparison demonstrates errors consistent with the project acceptance check.
If the set fails that check, the smallest evidence-based next action is to
obtain the matching 25 degC HPPC/pulse data and identify the same 1-RC topology;
generic or manually tuned replacement values should not be substituted.
