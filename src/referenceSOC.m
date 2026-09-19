function SOC_ref = referenceSOC(t, I_clean, SOC0, Q_As)
%REFERENCESOC Integrate clean current into the reference SOC trajectory.
%   SOC_REF = REFERENCESOC(T, I_CLEAN, SOC0, Q_AS) takes numeric time (s),
%   clean sign-normalised current (A, positive = discharge), initial SOC,
%   and capacity in amp-seconds.
%
%   Integrates every real interval using the right-endpoint convention and
%   ideal coulombic efficiency:
%     SOC_ref(1) = SOC0
%     SOC_ref(k) = SOC_ref(k-1) - I_clean(k)*(t(k)-t(k-1))/Q_As
%
%   No interpolation, no constant-dt assumption, no touching voltage or
%   dataset Ah, no clamping. If SOC drifts outside [0, 1] it stays that
%   way here - that's a downstream problem, not this function's job to
%   hide.

arguments
    t (:, 1) double
    I_clean (:, 1) double
    SOC0 (1, 1) double
    Q_As (1, 1) double
end

assert(~isempty(t), "referenceSOC:EmptyInput", ...
    "Time and current vectors must not be empty.");
assert(numel(t) == numel(I_clean), "referenceSOC:SizeMismatch", ...
    "Time and current vectors must contain the same number of samples.");
assert(all(isfinite(t)) && all(isfinite(I_clean)), ...
    "referenceSOC:NonFiniteVector", ...
    "Time and current vectors must contain only finite values.");
assert(all(diff(t) > 0), "referenceSOC:InvalidTime", ...
    "Numeric time must be strictly increasing.");
assert(isfinite(SOC0) && SOC0 >= 0 && SOC0 <= 1, ...
    "referenceSOC:InvalidInitialSOC", ...
    "Initial SOC must be a finite scalar in [0, 1].");
assert(isfinite(Q_As) && Q_As > 0, "referenceSOC:InvalidCapacity", ...
    "Capacity must be a positive finite scalar in ampere-seconds.");

SOC_ref = zeros(size(t));
SOC_ref(1) = SOC0;
if numel(t) > 1
    socChange = I_clean(2:end) .* diff(t) / Q_As;
    SOC_ref(2:end) = SOC0 - cumsum(socChange);
end
end
