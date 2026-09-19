function tests = testEstimateHPPC1RC
%TESTIMATEHPPC1RC Tests for the approved one-pass HPPC calculation.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
end

function testSyntheticFixed64WindowRecoveryAndProvenance(testCase)
trueR0_Ohm = 0.012;
trueR1_Ohm = 0.034;
trueTau_s = 3.7;
fitView = syntheticFitView(trueR0_Ohm, trueR1_Ohm, trueTau_s);

estimate = estimateHPPC1RC(fitView);

verifyEqual(testCase, estimate.R0, trueR0_Ohm, 'RelTol', 1e-5);
verifyEqual(testCase, estimate.R1, trueR1_Ohm, 'RelTol', 1e-5);
verifyEqual(testCase, estimate.tau, trueTau_s, 'RelTol', 1e-5);
verifyEqual(testCase, estimate.C1, estimate.tau / estimate.R1, 'AbsTol', 1e-14);
verifyTrue(testCase, all(isfinite([estimate.R0, estimate.R1, estimate.tau, estimate.C1])));
verifyTrue(testCase, all([estimate.R0, estimate.R1, estimate.tau, estimate.C1] > 0));
verifyGreaterThan(testCase, estimate.method.exitFlag, 0);
verifyEqual(testCase, estimate.method.candidateWindowCount, 64);
verifyEqual(testCase, numel(estimate.windows), 64);
verifyEqual(testCase, estimate.residual.sampleCount, 640);
verifyLessThan(testCase, estimate.residual.maxAbs_V, 1e-8);
verifyEqual(testCase, estimate.residual.definition, ...
    "predicted voltage drop - observed local voltage drop");

for windowIndex = 1:64
    window = estimate.windows(windowIndex);
    verifyEqual(testCase, window.pulseSourceRows, window.pulseFitViewIndices);
    verifyEqual(testCase, window.precedingRestSourceRows, window.precedingRestFitViewIndices);
    verifyEqual(testCase, window.followingRestSourceRows, window.followingRestFitViewIndices);
    verifyEqual(testCase, window.baselineSourceRow, window.baselineFitViewIndex);
    verifyEqual(testCase, window.pulseSampleCount, 10);
    verifyGreaterThanOrEqual(testCase, window.pulseDuration_s, 9);
end
end

function testOptimisationFailureIsRejected(testCase)
fitView = syntheticFitView(0.012, 0.034, 3.7);
fitView.V(fitView.I < 0.05) = realmax;
fitView.V(fitView.I >= 0.05) = -realmax;

verifyError(testCase, @() estimateHPPC1RC(fitView), ...
    "estimateHPPC1RC:OptimisationFailure");
end

function testRealHPPCReturnsApproved64WindowsAndResidualEvidence(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
hppcFile = fullfile(projectRoot, 'data', 'raw', ...
    '03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat');
loaded = load(hppcFile, 'meas');
fitView = prepareHPPCFitView(loaded.meas);

estimate = estimateHPPC1RC(fitView);

verifyEqual(testCase, numel(estimate.windows), 64);
verifyEqual(testCase, estimate.method.candidateWindowCount, 64);
verifyEqual(testCase, estimate.method.exitFlag, 1);
verifyTrue(testCase, all(isfinite([estimate.R0, estimate.R1, estimate.tau, estimate.C1])));
verifyTrue(testCase, all([estimate.R0, estimate.R1, estimate.tau, estimate.C1] > 0));
verifyEqual(testCase, estimate.C1, estimate.tau / estimate.R1, 'AbsTol', 1e-14);
verifyEqual(testCase, estimate.R0, 0.0312244650416, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.R1, 0.0113287933401, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.tau, 5.23962539484, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.C1, 462.505161632, 'RelTol', 1e-8);
verifyTrue(testCase, all(isfinite(estimate.residual.values_V)));
verifyEqual(testCase, estimate.residual.sampleCount, numel(estimate.residual.values_V));
verifyEqual(testCase, estimate.residual.sumSquared_V2, ...
    sum(estimate.residual.values_V .^ 2), 'AbsTol', 1e-12);
verifyEqual(testCase, estimate.residual.RMSE_V, ...
    sqrt(mean(estimate.residual.values_V .^ 2)), 'AbsTol', 1e-14);
verifyEqual(testCase, estimate.residual.RMSE_V, 0.0763801145593, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.residual.MAE_V, 0.0397436893111, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.residual.bias_V, -0.0138694662193, 'RelTol', 1e-8);
verifyEqual(testCase, estimate.residual.maxAbs_V, 0.409215005265, 'RelTol', 1e-8);

for windowIndex = 1:numel(estimate.windows)
    window = estimate.windows(windowIndex);
    verifyGreaterThanOrEqual(testCase, window.pulseDuration_s, 9);
    verifyEqual(testCase, window.baselineFitViewIndex, ...
        window.precedingRestFitViewIndices(end));
    verifyEqual(testCase, window.baselineSourceRow, ...
        window.precedingRestSourceRows(end));
    verifyTrue(testCase, all(fitView.I(window.pulseFitViewIndices) >= 0.05));
    verifyTrue(testCase, all(fitView.I(window.precedingRestFitViewIndices) < 0.05));
    verifyTrue(testCase, all(fitView.I(window.followingRestFitViewIndices) < 0.05));
    verifyEqual(testCase, window.pulseSourceRows, ...
        fitView.sourceRow(window.pulseFitViewIndices));
end
end

function fitView = syntheticFitView(R0_Ohm, R1_Ohm, tau_s)
windowCount = 64;
samplesPerPulse = 10;
samplesPerWindow = samplesPerPulse + 2;
sampleCount = windowCount * samplesPerWindow;

time_s = zeros(sampleCount, 1);
current_A = zeros(sampleCount, 1);
voltage_V = zeros(sampleCount, 1);
for windowIndex = 1:windowCount
    firstIndex = (windowIndex - 1) * samplesPerWindow + 1;
    restBeforeIndex = firstIndex;
    pulseIndices = (firstIndex + 1:firstIndex + samplesPerPulse).';
    restAfterIndex = firstIndex + samplesPerWindow - 1;
    localTime_s = (0:samplesPerWindow - 1).';
    time_s(firstIndex:restAfterIndex) = localTime_s + (windowIndex - 1) * 20;
    pulseCurrent_A = 0.8 + 0.1 * mod(windowIndex - 1, 5);
    current_A(pulseIndices) = pulseCurrent_A;
    baselineVoltage_V = 4.10 - 0.005 * windowIndex;
    voltage_V(restBeforeIndex) = baselineVoltage_V;

    polarisation_V = 0;
    previousIndex = restBeforeIndex;
    for sampleIndex = pulseIndices.'
        dt_s = time_s(sampleIndex) - time_s(previousIndex);
        alpha = exp(-dt_s / tau_s);
        polarisation_V = alpha * polarisation_V + ...
            (1 - alpha) * R1_Ohm * current_A(sampleIndex);
        voltage_V(sampleIndex) = baselineVoltage_V - ...
            R0_Ohm * current_A(sampleIndex) - polarisation_V;
        previousIndex = sampleIndex;
    end
    voltage_V(restAfterIndex) = baselineVoltage_V;
end

fitView = struct;
fitView.t = time_s;
fitView.I = current_A;
fitView.V = voltage_V;
fitView.sourceRow = (1:sampleCount).';
end
