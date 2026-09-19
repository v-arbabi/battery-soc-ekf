function [V_model, V1] = ecmSimulate(t, I, SOC_ref, evaluateOCV, params)
%ECMSIMULATE Simulate 1-RC ECM voltage on the right-endpoint convention.
%   [V_MODEL, V1] = ECMSIMULATE(T, I, SOC_REF, EVALUATEOCV, PARAMS) takes
%   numeric time (s), sign-normalised current (A, positive = discharge),
%   and an externally supplied reference-SOC trajectory. EVALUATEOCV
%   evaluates OCV (V) at a given SOC. PARAMS needs R0, R1, tau.
%
%   V1(1) starts at zero polarisation. V_model(1) is NaN - there's no
%   preceding interval for the first sample, so nothing to compute yet.
%   For k >= 2:
%     dt_k = t(k) - t(k-1)
%     alpha_k = exp(-dt_k / tau)
%     V1(k) = alpha_k*V1(k-1) + (1-alpha_k)*R1*I(k)
%     V_model(k) = OCV(SOC_ref(k)) - R0*I(k) - V1(k)
%
%   Purely a forward simulator: it takes SOC_ref as given, doesn't fit or
%   tune anything, and never looks at the measured voltage.

arguments
    t (:, 1) double
    I (:, 1) double
    SOC_ref (:, 1) double
    evaluateOCV (1, 1) function_handle
    params (1, 1) struct
end

assert(~isempty(t), "ecmSimulate:EmptyInput", ...
    "Time, current, and reference-SOC vectors must not be empty.");
assert(numel(t) == numel(I) && numel(t) == numel(SOC_ref), ...
    "ecmSimulate:SizeMismatch", ...
    "Time, current, and reference-SOC vectors must have equal lengths.");
assert(all(isfinite(t)) && all(isfinite(I)) && all(isfinite(SOC_ref)), ...
    "ecmSimulate:NonFiniteVector", ...
    "Time, current, and reference-SOC vectors must contain finite values.");
assert(all(SOC_ref >= 0 & SOC_ref <= 1), ...
    "ecmSimulate:InvalidReferenceSOC", ...
    "Reference SOC must remain within the approved physical interval [0, 1].");
assert(all(diff(t) > 0), "ecmSimulate:InvalidTime", ...
    "Numeric time must be strictly increasing.");
assert(all(isfield(params, {'R0', 'R1', 'tau'})), ...
    "ecmSimulate:MissingParameter", ...
    "Parameters must provide R0, R1, and tau.");
assert(isscalar(params.R0) && isscalar(params.R1) && isscalar(params.tau), ...
    "ecmSimulate:InvalidParameter", ...
    "R0, R1, and tau must be scalar values.");

parameterValues = [params.R0, params.R1, params.tau];
assert(isnumeric(parameterValues) && isreal(parameterValues) && ...
    all(isfinite(parameterValues)) && all(parameterValues > 0), ...
    "ecmSimulate:InvalidParameter", ...
    "R0, R1, and tau must be positive finite real values.");

V1 = zeros(size(t));
V_model = nan(size(t));

for k = 2:numel(t)
    dt_s = t(k) - t(k - 1);
    alpha = exp(-dt_s / params.tau);
    V1(k) = alpha * V1(k - 1) + (1 - alpha) * params.R1 * I(k);
end

if numel(t) > 1
    ocv_V = evaluateOCV(SOC_ref(2:end));
    assert(isnumeric(ocv_V) && isreal(ocv_V) && ...
        isequal(size(ocv_V), size(SOC_ref(2:end))) && all(isfinite(ocv_V)), ...
        "ecmSimulate:InvalidOCVOutput", ...
        "The OCV evaluator must return a finite real vector matching its input.");
    V_model(2:end) = ocv_V - params.R0 * I(2:end) - V1(2:end);
end
end
