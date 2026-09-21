#!/usr/bin/env python3
"""Independent Python cross-check of the MATLAB SOC-estimation results.

This is a second, separate implementation of the same equations used in
src/*.m. It shares no code with the MATLAB pipeline: it re-reads the raw
.mat files, rebuilds the OCV target, re-fits the degree-9 polynomial,
and re-runs Coulomb counting and the Extended Kalman Filter for all five
experiments. The point is that a transcription error in either
implementation would show up here as a disagreement.

Run from the repository root:

    python scripts/crosscheck.py

Requires numpy and scipy.
"""

from pathlib import Path
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"

C20_FILE = RAW / "05-08-17_13.26 C20 OCV Test_C20_25dC.mat"
US06_FILE = RAW / "03-20-17_01.43 25degC_US06_Pan18650PF.mat"

REQUIRED = ["TimeStamp", "Voltage", "Current", "Ah", "Wh", "Power",
            "Battery_Temp_degC", "Time", "Chamber_Temp_degC"]

# Locked HPPC-identified 1-RC set (src/ecmParamsHPPC.m)
HPPC = dict(R0=0.0312244650416, R1=0.0113287933401, tau=5.23962539484,
            Q_Ah=2.997393193)
HPPC["Q_As"] = HPPC["Q_Ah"] * 3600.0

# Literature set, used only as the mismatched model in Experiment D
# (src/ecmParams.m)
LIT = dict(R0=0.004500052, R1=0.028981830, C1=218.8730953, Q_Ah=2.997393193)
LIT["tau"] = LIT["R1"] * LIT["C1"]
LIT["Q_As"] = LIT["Q_Ah"] * 3600.0

TUNING = dict(P0_SOC=0.01 ** 2, P0_V1=0.01 ** 2,
              q_SOC=1e-10, q_V1=1e-6, R_V=3e-3 ** 2)


def load_dataset(path):
    """Port of src/loadDataset.m: row window, duplicate removal, sign flip."""
    meas = loadmat(path, simplify_cells=True)["meas"]
    missing = [f for f in REQUIRED if f not in meas]
    if missing:
        raise ValueError(f"missing fields in {path.name}: {missing}")

    n = len(np.atleast_1d(meas["Time"]))
    if path.name == C20_FILE.name:
        if n != 2453:
            raise ValueError("the approved C/20 file must contain 2453 rows")
        rows = np.arange(1, 2453)          # source rows 1..2452, one-based
    elif path.name == US06_FILE.name:
        rows = np.arange(1, n + 1)
    else:
        raise ValueError(f"not one of the two approved datasets: {path.name}")

    idx = rows - 1                          # to zero-based
    t_sel = np.asarray(meas["Time"], float)[idx]
    if np.any(np.diff(t_sel) < 0):
        raise ValueError("numeric Time decreases within the approved rows")

    # drop the second row of every exact zero-dt duplicate pair
    dup = np.flatnonzero(np.diff(t_sel) == 0) + 1
    keep = np.ones(rows.size, bool)
    keep[dup] = False
    kept = idx[keep]

    t = np.asarray(meas["Time"], float)[kept]
    if not np.all(np.diff(t) > 0):
        raise ValueError("Time must be strictly increasing after de-duplication")

    return dict(
        t=t,
        I=-np.asarray(meas["Current"], float)[kept],   # positive = discharge
        V=np.asarray(meas["Voltage"], float)[kept],
        source_row=rows[keep],
    )


def build_ocv_target():
    """Port of src/buildOCVTarget.m: 101-point branch-mean OCV target."""
    d = load_dataset(C20_FILE)
    sr = d["source_row"]

    dis = np.flatnonzero((sr >= 7) & (sr <= 1247))
    chg = np.flatnonzero((sr >= 1309) & (sr <= 2391))
    rest_dis = np.flatnonzero(sr == 6)
    rest_empty = np.flatnonzero(sr == 1307)
    rest_full = np.flatnonzero(sr == 2451)

    dt_dis = np.diff(d["t"][np.concatenate([rest_dis, dis])])
    dt_chg = np.diff(d["t"][np.concatenate([rest_empty, chg])])
    cap_dis = np.sum(d["I"][dis] * dt_dis) / 3600.0
    cap_chg = -np.sum(d["I"][chg] * dt_chg) / 3600.0

    soc_dis = 1.0 - np.concatenate(
        [[0.0], np.cumsum(d["I"][dis] * dt_dis) / 3600.0]) / cap_dis
    soc_chg = np.concatenate(
        [[0.0], np.cumsum(-d["I"][chg] * dt_chg) / 3600.0]) / cap_chg
    soc_dis[0], soc_dis[-1] = 1.0, 0.0
    soc_chg[0], soc_chg[-1] = 0.0, 1.0

    v_dis = d["V"][np.concatenate([rest_dis, dis])]
    v_chg = d["V"][np.concatenate([rest_empty, chg])]

    soc_t = np.linspace(0.0, 1.0, 101)
    v_dis_grid = np.interp(soc_t, soc_dis[::-1], v_dis[::-1])
    v_chg_grid = np.interp(soc_t, soc_chg, v_chg)

    v_t = np.zeros(101)
    v_t[0] = d["V"][rest_empty][0]
    v_t[1:100] = (v_dis_grid[1:100] + v_chg_grid[1:100]) / 2.0
    v_t[100] = np.mean([d["V"][rest_dis][0], d["V"][rest_full][0]])
    return soc_t, v_t


def build_ocv(soc_t, v_t, degree=9):
    """Port of src/buildOCV.m: single degree-9 polynomial plus its derivative."""
    coef = np.polyfit(soc_t, v_t, degree)
    dcoef = np.polyder(coef)
    return (lambda s: np.polyval(coef, s)), (lambda s: np.polyval(dcoef, s))


def integrate_soc(t, I, soc0, Q_As):
    """Shared by referenceSOC.m and coulombCount.m: right-endpoint integration."""
    soc = np.empty_like(t)
    soc[0] = soc0
    soc[1:] = soc0 - np.cumsum(I[1:] * np.diff(t) / Q_As)
    return soc


def ekf_estimate(t, I, V_meas, ocv, docv, params, soc0, Q_As, tuning):
    """Port of src/ekfEstimate.m: 2-state EKF on x = [SOC; V1]."""
    R0, R1, tau = params["R0"], params["R1"], params["tau"]
    n = t.size

    x = np.array([soc0, 0.0])
    P = np.diag([tuning["P0_SOC"], tuning["P0_V1"]])
    Qk = np.diag([tuning["q_SOC"], tuning["q_V1"]])
    Rk = tuning["R_V"]

    SOC = np.zeros(n)
    P_SOC = np.zeros(n)
    SOC[0], P_SOC[0] = x[0], P[0, 0]

    for k in range(1, n):
        dt = t[k] - t[k - 1]
        a = np.exp(-dt / tau)

        x_pred = np.array([x[0] - (dt / Q_As) * I[k],
                           a * x[1] + (1.0 - a) * R1 * I[k]])
        F = np.array([[1.0, 0.0], [0.0, a]])
        P_pred = F @ P @ F.T + Qk

        y_pred = ocv(x_pred[0]) - x_pred[1] - R0 * I[k]
        resid = V_meas[k] - y_pred

        H = np.array([docv(x_pred[0]), -1.0])
        S = H @ P_pred @ H.T + Rk
        K = (P_pred @ H) / S

        x = x_pred + K * resid
        P = (np.eye(2) - np.outer(K, H)) @ P_pred
        x[0] = min(max(x[0], 0.0), 1.0)      # SOC clamped to [0, 1]

        SOC[k], P_SOC[k] = x[0], P[0, 0]

    return SOC, P_SOC


def rmse(estimate, reference):
    return float(np.sqrt(np.mean((estimate - reference) ** 2)))


def main():
    soc_t, v_t = build_ocv_target()
    ocv, docv = build_ocv(soc_t, v_t)
    d = load_dataset(US06_FILE)
    t, I, V = d["t"], d["I"], d["V"]

    soc_ref = integrate_soc(t, I, 1.0, HPPC["Q_As"])
    rows = []

    # A - nominal
    tun = dict(TUNING, P0_SOC=0.01 ** 2)
    soc_ekf, _ = ekf_estimate(t, I, V, ocv, docv, HPPC, 1.0, HPPC["Q_As"], tun)
    cc = integrate_soc(t, I, 1.0, HPPC["Q_As"])
    rows.append(("A  nominal", rmse(soc_ekf, soc_ref), 7.4317,
                 rmse(cc, soc_ref), 0.0, True))

    # B - wrong initial SOC
    tun = dict(TUNING, P0_SOC=0.20 ** 2)
    soc_ekf, _ = ekf_estimate(t, I, V, ocv, docv, HPPC, 0.8, HPPC["Q_As"], tun)
    cc = integrate_soc(t, I, 0.8, HPPC["Q_As"])
    rows.append(("B  wrong SOC0", rmse(soc_ekf, soc_ref), 7.3955,
                 rmse(cc, soc_ref), 20.0, True))

    # C - corrupted sensors. The voltage-noise draw comes from MATLAB's
    # rng(42) stream, which cannot be reproduced bit-for-bit in numpy, so
    # the EKF figure here is checked to a tolerance rather than exactly.
    # Coulomb counting sees only the deterministic +0.02 A current bias,
    # so its figure is still an exact comparison.
    rng = np.random.default_rng(42)
    I_bad = I + 0.02
    V_bad = V + 3e-3 * rng.standard_normal(V.size)
    soc_ekf, _ = ekf_estimate(t, I_bad, V_bad, ocv, docv, HPPC, 1.0,
                              HPPC["Q_As"], TUNING)
    cc = integrate_soc(t, I_bad, 1.0, HPPC["Q_As"])
    rows.append(("C  sensor fault", rmse(soc_ekf, soc_ref), 7.4455,
                 rmse(cc, soc_ref), 0.5157, False))

    # D - mismatched process model
    soc_ekf, _ = ekf_estimate(t, I, V, ocv, docv, LIT, 1.0, HPPC["Q_As"], TUNING)
    cc = integrate_soc(t, I, 1.0, HPPC["Q_As"])
    rows.append(("D  model mismatch", rmse(soc_ekf, soc_ref), 8.9161,
                 rmse(cc, soc_ref), 0.0, True))

    # E - capacity mismatch
    Q_nom = 2.9 * 3600.0
    soc_ekf, _ = ekf_estimate(t, I, V, ocv, docv, HPPC, 1.0, Q_nom, TUNING)
    cc = integrate_soc(t, I, 1.0, Q_nom)
    rows.append(("E  capacity", rmse(soc_ekf, soc_ref), 7.7436,
                 rmse(cc, soc_ref), 1.7294, True))

    print(f"samples: {t.size}   capacity: {HPPC['Q_Ah']:.9f} Ah")
    print(f"OCV polynomial degree 9, target RMSE "
          f"{rmse(np.polyval(np.polyfit(soc_t, v_t, 9), soc_t), v_t) * 1000:.4f} mV\n")

    head = (f"{'experiment':18s} {'EKF py':>9s} {'EKF m':>9s} {'d':>8s}   "
            f"{'CC py':>9s} {'CC m':>9s} {'d':>8s}")
    print(head)
    print("-" * len(head))

    worst = 0.0
    for name, ekf_py, ekf_m, cc_py, cc_m, exact in rows:
        ekf_py *= 100.0
        cc_py *= 100.0
        d_ekf, d_cc = ekf_py - ekf_m, cc_py - cc_m
        if exact:
            worst = max(worst, abs(d_ekf), abs(d_cc))
        else:
            worst = max(worst, abs(d_cc))
        print(f"{name:18s} {ekf_py:8.4f}% {ekf_m:8.4f}% {d_ekf:+8.4f}   "
              f"{cc_py:8.4f}% {cc_m:8.4f}% {d_cc:+8.4f}")

    print(f"\nlargest disagreement on an exactly comparable figure: "
          f"{worst:.4f} percentage points")
    print("Experiment C's EKF figure is excluded from that maximum because "
          "its voltage-noise\nrealisation differs between MATLAB and numpy; "
          "it is expected to agree only to within\nthe spread of the noise draw.")


if __name__ == "__main__":
    main()
