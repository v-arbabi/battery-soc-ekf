function SOC_CC = coulombCount(t, I, SOC0, Q_As)
%COULOMBCOUNT Integrate current into the SOC_CC baseline estimate.
%   SOC_CC = COULOMBCOUNT(T, I, SOC0, Q_AS) takes numeric time (s),
%   sign-normalised current (A, positive = discharge), the initial SOC,
%   and an assumed capacity in amp-seconds.
%
%   Uses every real interval in the data (no resampling) with ideal
%   coulombic efficiency and the right-endpoint convention:
%     SOC_CC(1) = SOC0
%     SOC_CC(k) = SOC_CC(k-1) - I(k)*(t(k)-t(k-1))/Q_As
%
%   Kept separate from SOC_ref and the EKF on purpose: this is the raw,
%   uncorrected coulomb-counting baseline, so it doesn't touch voltage,
%   OCV, or dataset Ah, and it doesn't clamp to [0,1]. If it drifts
%   outside that range that's useful diagnostic info, so it's left alone.

arguments
    t (:, 1) double
    I (:, 1) double
    SOC0 (1, 1) double
    Q_As (1, 1) double
end

assert(~isempty(t), "coulombCount:EmptyInput", ...
    "Time and current vectors must not be empty.");
assert(numel(t) == numel(I), "coulombCount:SizeMismatch", ...
    "Time and current vectors must contain the same number of samples.");
assert(all(isfinite(t)) && all(isfinite(I)), ...
    "coulombCount:NonFiniteVector", ...
    "Time and current vectors must contain only finite values.");
assert(all(diff(t) > 0), "coulombCount:InvalidTime", ...
    "Numeric time must be strictly increasing.");
assert(isfinite(SOC0) && SOC0 >= 0 && SOC0 <= 1, ...
    "coulombCount:InvalidInitialSOC", ...
    "Initial SOC must be a finite scalar in [0, 1].");
assert(isfinite(Q_As) && Q_As > 0, "coulombCount:InvalidCapacity", ...
    "Capacity must be a positive finite scalar in ampere-seconds.");

SOC_CC = zeros(size(t));
SOC_CC(1) = SOC0;
if numel(t) > 1
    socChange = I(2:end) .* diff(t) / Q_As;
    SOC_CC(2:end) = SOC0 - cumsum(socChange);
end
end
