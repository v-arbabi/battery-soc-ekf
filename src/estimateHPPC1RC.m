function estimate = estimateHPPC1RC(fitView)
%ESTIMATEHPPC1RC Run the one-pass HPPC 1-RC parameter calculation.
%   ESTIMATE = ESTIMATEHPPC1RC(FITVIEW) works only on the strict-time HPPC
%   fit view. It picks up every active pulse lasting at least 9 s, the
%   fixed 64-window set approved for this calculation, and for each one
%   the last sample of the preceding rest becomes the local voltage
%   baseline, with the pulse-local RC drop starting from zero.
%
%   Right-endpoint pulse recurrence (positive = discharge):
%     V1(k) = exp(-dt/tau)*V1(k-1) + (1-exp(-dt/tau))*R1*I(k)
%     predictedDrop(k) = R0*I(k) + V1(k)
%   The first pulse sample uses the real interval back to the preceding
%   rest endpoint. Parameters come from one global least-squares fit in
%   log(R0), log(R1), log(tau) (keeps everything positive by
%   construction), with C1 = tau/R1 derived afterward. Doesn't touch
%   ecmParams.m and doesn't adopt these values automatically - that's a
%   separate sign-off.

arguments
    fitView (1, 1) struct
end

requiredFields = ["t", "I", "V", "sourceRow"];
assert(all(isfield(fitView, requiredFields)), ...
    "estimateHPPC1RC:MissingField", ...
    "The fit view must provide t, I, V, and sourceRow.");

time_s = fitView.t(:);
current_A = fitView.I(:);
voltage_V = fitView.V(:);
sourceRow = fitView.sourceRow(:);
sampleCount = numel(time_s);

assert(sampleCount > 1, "estimateHPPC1RC:InsufficientSamples", ...
    "The fit view must contain at least two retained samples.");
assert(numel(current_A) == sampleCount && numel(voltage_V) == sampleCount && ...
    numel(sourceRow) == sampleCount, "estimateHPPC1RC:SizeMismatch", ...
    "Fit-view time, current, voltage, and source rows must have equal lengths.");
assert(all(isfinite(time_s)) && all(isfinite(current_A)) && all(isfinite(voltage_V)), ...
    "estimateHPPC1RC:NonFiniteSignal", ...
    "Fit-view time, current, and voltage must be finite.");
assert(all(diff(time_s) > 0), "estimateHPPC1RC:InvalidTime", ...
    "Fit-view time must be strictly increasing.");
assert(all(sourceRow == floor(sourceRow)) && all(sourceRow > 0), ...
    "estimateHPPC1RC:InvalidSourceRows", ...
    "Fit-view source rows must be positive integer indices.");

activeThreshold_A = 0.05;
fullPulseCoverageThreshold_s = 9.0;
isActive = current_A >= activeThreshold_A;
segmentStart = [1; find(diff(isActive) ~= 0) + 1];
segmentEnd = [segmentStart(2:end) - 1; sampleCount];
segmentIsActive = isActive(segmentStart);
pulseSegmentIndices = find(segmentIsActive);

windows = repmat(localEmptyWindow(), 0, 1);
for segmentIndex = pulseSegmentIndices.'
    pulseIndices = (segmentStart(segmentIndex):segmentEnd(segmentIndex)).';
    pulseDuration_s = time_s(pulseIndices(end)) - time_s(pulseIndices(1));
    hasPrecedingRest = segmentIndex > 1 && ~segmentIsActive(segmentIndex - 1);
    hasFollowingRest = segmentIndex < numel(segmentIsActive) && ...
        ~segmentIsActive(segmentIndex + 1);

    if ~(hasPrecedingRest && hasFollowingRest && ...
            pulseDuration_s >= fullPulseCoverageThreshold_s)
        continue;
    end

    precedingRestIndices = ...
        (segmentStart(segmentIndex - 1):segmentEnd(segmentIndex - 1)).';
    followingRestIndices = ...
        (segmentStart(segmentIndex + 1):segmentEnd(segmentIndex + 1)).';
    assert(all(current_A(pulseIndices) >= activeThreshold_A) && ...
        all(current_A(precedingRestIndices) < activeThreshold_A) && ...
        all(current_A(followingRestIndices) < activeThreshold_A), ...
        "estimateHPPC1RC:InvalidWindow", ...
        "The fixed candidate segmentation must retain pulse and rest states.");

    window = localEmptyWindow();
    window.pulseFitViewIndices = pulseIndices;
    window.pulseSourceRows = sourceRow(pulseIndices);
    window.precedingRestFitViewIndices = precedingRestIndices;
    window.precedingRestSourceRows = sourceRow(precedingRestIndices);
    window.followingRestFitViewIndices = followingRestIndices;
    window.followingRestSourceRows = sourceRow(followingRestIndices);
    window.baselineFitViewIndex = precedingRestIndices(end);
    window.baselineSourceRow = sourceRow(precedingRestIndices(end));
    window.baselineVoltage_V = voltage_V(precedingRestIndices(end));
    window.pulseStartTime_s = time_s(pulseIndices(1));
    window.pulseEndTime_s = time_s(pulseIndices(end));
    window.pulseDuration_s = pulseDuration_s;
    window.pulseSampleCount = numel(pulseIndices);
    windows(end + 1, 1) = window; %#ok<AGROW>
end

assert(numel(windows) == 64, "estimateHPPC1RC:UnexpectedWindowCount", ...
    "The approved calculation requires exactly 64 fixed full-coverage pulses.");

initialParameters = struct('R0', 0.01, 'R1', 0.01, 'tau', 10.0);
initialLogParameters = log([initialParameters.R0; initialParameters.R1; initialParameters.tau]);
objective = @(logParameters) localSumSquaredResiduals( ...
    logParameters, windows, time_s, current_A, voltage_V);
optimisationOptions = optimset('Display', 'off', 'MaxIter', 2000, ...
    'MaxFunEvals', 6000, 'TolX', 1e-10, 'TolFun', 1e-12);
try
    [logParameters, sumSquaredResidual_V2, exitFlag, optimisationOutput] = ...
        fminsearch(objective, initialLogParameters, optimisationOptions);
catch optimisationException
    error("estimateHPPC1RC:OptimisationFailure", ...
        "The fixed global optimisation failed: %s", optimisationException.message);
end
assert(exitFlag > 0 && isfinite(sumSquaredResidual_V2) && ...
    all(isfinite(logParameters)), "estimateHPPC1RC:OptimisationFailure", ...
    "The fixed global optimisation did not terminate successfully.");

parameterValues = exp(logParameters);
R0_Ohm = parameterValues(1);
R1_Ohm = parameterValues(2);
tau_s = parameterValues(3);
C1_F = tau_s / R1_Ohm;
candidateParameterValues = [R0_Ohm, R1_Ohm, tau_s, C1_F];
assert(isreal(candidateParameterValues) && all(isfinite(candidateParameterValues)) && ...
    all(candidateParameterValues > 0), "estimateHPPC1RC:OptimisationFailure", ...
    "Transformed candidate parameters must be finite and strictly positive.");

[residual_V, predictedDrop_V, observedDrop_V, sampleWindowIndex] = ...
    localResiduals([R0_Ohm; R1_Ohm; tau_s], windows, time_s, current_A, voltage_V);

for windowIndex = 1:numel(windows)
    isWindowSample = sampleWindowIndex == windowIndex;
    windowResidual_V = residual_V(isWindowSample);
    windows(windowIndex).predictedDrop_V = predictedDrop_V(isWindowSample);
    windows(windowIndex).observedDrop_V = observedDrop_V(isWindowSample);
    windows(windowIndex).residual_V = windowResidual_V;
    windows(windowIndex).residualRMSE_V = sqrt(mean(windowResidual_V .^ 2));
    windows(windowIndex).residualBias_V = mean(windowResidual_V);
    windows(windowIndex).residualMaxAbs_V = max(abs(windowResidual_V));
end

estimate = struct;
estimate.R0 = R0_Ohm;
estimate.R1 = R1_Ohm;
estimate.tau = tau_s;
estimate.C1 = C1_F;
estimate.windows = windows;
estimate.residual = struct;
estimate.residual.definition = "predicted voltage drop - observed local voltage drop";
estimate.residual.sampleCount = numel(residual_V);
estimate.residual.sumSquared_V2 = sumSquaredResidual_V2;
estimate.residual.RMSE_V = sqrt(mean(residual_V .^ 2));
estimate.residual.MAE_V = mean(abs(residual_V));
estimate.residual.bias_V = mean(residual_V);
estimate.residual.maxAbs_V = max(abs(residual_V));
estimate.residual.values_V = residual_V;
estimate.method = struct;
estimate.method.activeThreshold_A = activeThreshold_A;
estimate.method.fullPulseCoverageThreshold_s = fullPulseCoverageThreshold_s;
estimate.method.candidateWindowCount = numel(windows);
estimate.method.baselinePolicy = "preceding retained-rest endpoint voltage";
estimate.method.initialRCState_V = 0;
estimate.method.timingConvention = "right endpoint using I(k) over [t(k-1), t(k)]";
estimate.method.optimisation = "one global fminsearch least-squares pass in log(R0), log(R1), log(tau)";
estimate.method.initialParameters = initialParameters;
estimate.method.exitFlag = exitFlag;
estimate.method.output = optimisationOutput;
end

function window = localEmptyWindow()
window = struct( ...
    'pulseFitViewIndices', zeros(0, 1), ...
    'pulseSourceRows', zeros(0, 1), ...
    'precedingRestFitViewIndices', zeros(0, 1), ...
    'precedingRestSourceRows', zeros(0, 1), ...
    'followingRestFitViewIndices', zeros(0, 1), ...
    'followingRestSourceRows', zeros(0, 1), ...
    'baselineFitViewIndex', NaN, ...
    'baselineSourceRow', NaN, ...
    'baselineVoltage_V', NaN, ...
    'pulseStartTime_s', NaN, ...
    'pulseEndTime_s', NaN, ...
    'pulseDuration_s', NaN, ...
    'pulseSampleCount', NaN, ...
    'predictedDrop_V', zeros(0, 1), ...
    'observedDrop_V', zeros(0, 1), ...
    'residual_V', zeros(0, 1), ...
    'residualRMSE_V', NaN, ...
    'residualBias_V', NaN, ...
    'residualMaxAbs_V', NaN);
end

function sumSquared_V2 = localSumSquaredResiduals(logParameters, windows, time_s, current_A, voltage_V)
parameters = exp(logParameters);
[residual_V, ~, ~, ~] = localResiduals(parameters, windows, time_s, current_A, voltage_V);
sumSquared_V2 = sum(residual_V .^ 2);
end

function [residual_V, predictedDrop_V, observedDrop_V, sampleWindowIndex] = ...
    localResiduals(parameters, windows, time_s, current_A, voltage_V)
R0_Ohm = parameters(1);
R1_Ohm = parameters(2);
tau_s = parameters(3);

residual_V = zeros(0, 1);
predictedDrop_V = zeros(0, 1);
observedDrop_V = zeros(0, 1);
sampleWindowIndex = zeros(0, 1);
for windowIndex = 1:numel(windows)
    pulseIndices = windows(windowIndex).pulseFitViewIndices;
    previousIndex = windows(windowIndex).baselineFitViewIndex;
    polarisation_V = 0;
    for sampleIndex = pulseIndices.'
        dt_s = time_s(sampleIndex) - time_s(previousIndex);
        alpha = exp(-dt_s / tau_s);
        polarisation_V = alpha * polarisation_V + ...
            (1 - alpha) * R1_Ohm * current_A(sampleIndex);
        predictedDrop_V(end + 1, 1) = ...
            R0_Ohm * current_A(sampleIndex) + polarisation_V; %#ok<AGROW>
        observedDrop_V(end + 1, 1) = ...
            windows(windowIndex).baselineVoltage_V - voltage_V(sampleIndex); %#ok<AGROW>
        sampleWindowIndex(end + 1, 1) = windowIndex; %#ok<AGROW>
        previousIndex = sampleIndex;
    end
end
residual_V = predictedDrop_V - observedDrop_V;
end
