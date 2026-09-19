function result = ekfEstimate(t, I, V_meas, evaluateOCV, evaluateDOCV, params, SOC0, Q_As, tuning)
%EKFESTIMATE Extended Kalman Filter SOC estimate for the 1-RC ECM.
%   RESULT = EKFESTIMATE(T, I, V_MEAS, EVALUATEOCV, EVALUATEDOCV, PARAMS,
%   SOC0, Q_AS, TUNING) estimates the two-state vector x = [SOC; V1] from
%   current and measured terminal voltage. Uses the same right-endpoint
%   timing convention as ecmSimulate.m and referenceSOC.m: for k >= 2 both
%   the prediction and the measurement at sample k use the interval
%   [t(k-1), t(k)] and I(k).
%
%   State-space model:
%     Prediction (process model):
%       SOC(k)  = SOC(k-1) - (dt_k/Q_As) * I(k)
%       V1(k)   = alpha_k*V1(k-1) + (1-alpha_k)*R1*I(k),  alpha_k = exp(-dt_k/tau)
%     Measurement model:
%       V(k) = OCV(SOC(k)) - V1(k) - R0*I(k)
%   State transition is linear (F is constant). The measurement model is
%   nonlinear through OCV(SOC), so H gets re-evaluated every step at the
%   predicted SOC:
%       F = [1, 0; 0, alpha_k]
%       H = [dOCV/dSOC(SOC_pred), -1]
%
%   TUNING fields (all scalar):
%     P0_SOC, P0_V1   initial state-covariance diagonal
%     q_SOC, q_V1     process-noise covariance diagonal (Q)
%     R_V             measurement-noise variance (R), volts^2
%
%   RESULT fields (column vectors, length numel(t), except sigma3):
%     SOC, V1         filtered state estimates
%     P_SOC           filtered SOC-variance history, P(1,1)
%     sigma3_SOC      3*sqrt(P_SOC), the plotted uncertainty half-width
%     innovation      V_meas(k) - predicted V(k); innovation(1) = NaN
%     coverage3Sigma  percent of samples within the 3-sigma band; needs
%                     ground truth (SOC_ref), which this function never
%                     sees, so it's left for the caller to compute
%
%   Never touches SOC_ref, dataset Ah, or anything about corrupted-input
%   handling, just current, voltage, and the parameters it's handed, same
%   as ecmSimulate.m and coulombCount.m. That's what keeps the test fair:
%   the filter can't cheat by peeking at the answer.

arguments
    t (:, 1) double
    I (:, 1) double
    V_meas (:, 1) double
    evaluateOCV (1, 1) function_handle
    evaluateDOCV (1, 1) function_handle
    params (1, 1) struct
    SOC0 (1, 1) double
    Q_As (1, 1) double
    tuning (1, 1) struct
end

assert(~isempty(t), "ekfEstimate:EmptyInput", ...
    "Time, current, and voltage vectors must not be empty.");
assert(numel(t) == numel(I) && numel(t) == numel(V_meas), ...
    "ekfEstimate:SizeMismatch", ...
    "Time, current, and voltage vectors must have equal lengths.");
assert(all(isfinite(t)) && all(isfinite(I)) && all(isfinite(V_meas)), ...
    "ekfEstimate:NonFiniteVector", ...
    "Time, current, and voltage vectors must contain finite values.");
assert(all(diff(t) > 0), "ekfEstimate:InvalidTime", ...
    "Numeric time must be strictly increasing.");
assert(isfinite(SOC0) && SOC0 >= 0 && SOC0 <= 1, ...
    "ekfEstimate:InvalidInitialSOC", ...
    "Initial SOC must be a finite scalar in [0, 1].");
assert(isfinite(Q_As) && Q_As > 0, "ekfEstimate:InvalidCapacity", ...
    "Capacity must be a positive finite scalar in ampere-seconds.");
assert(all(isfield(params, {'R0', 'R1', 'tau'})), ...
    "ekfEstimate:MissingParameter", "Parameters must provide R0, R1, and tau.");
requiredTuning = ["P0_SOC", "P0_V1", "q_SOC", "q_V1", "R_V"];
assert(all(isfield(tuning, requiredTuning)), ...
    "ekfEstimate:MissingTuning", ...
    "Tuning must provide P0_SOC, P0_V1, q_SOC, q_V1, and R_V.");

N = numel(t);
R0 = params.R0;
R1 = params.R1;
tau = params.tau;

x = [SOC0; 0];
P = diag([tuning.P0_SOC, tuning.P0_V1]);
Qk = diag([tuning.q_SOC, tuning.q_V1]);
Rk = tuning.R_V;

SOC = zeros(N, 1);
V1 = zeros(N, 1);
P_SOC = zeros(N, 1);
innovation = nan(N, 1);

SOC(1) = x(1);
V1(1) = x(2);
P_SOC(1) = P(1, 1);

for k = 2:N
    dt_k = t(k) - t(k - 1);
    alpha = exp(-dt_k / tau);

    % predict
    x_pred = [ ...
        x(1) - (dt_k / Q_As) * I(k); ...
        alpha * x(2) + (1 - alpha) * R1 * I(k) ...
    ];
    F = [1, 0; 0, alpha];
    P_pred = F * P * F' + Qk;

    % update
    ocv_V = evaluateOCV(x_pred(1));
    y_pred = ocv_V - x_pred(2) - R0 * I(k);
    resid = V_meas(k) - y_pred;

    H = [evaluateDOCV(x_pred(1)), -1];
    S = H * P_pred * H' + Rk;
    K = (P_pred * H') / S;

    x_upd = x_pred + K * resid;
    P_upd = (eye(2) - K * H) * P_pred;

    % clamp to [0,1] - SOC can't go negative or over 100%
    x_upd(1) = min(max(x_upd(1), 0), 1);

    x = x_upd;
    P = P_upd;

    SOC(k) = x(1);
    V1(k) = x(2);
    P_SOC(k) = P(1, 1);
    innovation(k) = resid;
end

result = struct;
result.SOC = SOC;
result.V1 = V1;
result.P_SOC = P_SOC;
result.sigma3_SOC = 3 * sqrt(P_SOC);
result.innovation = innovation;
result.params = params;
result.tuning = tuning;
end
