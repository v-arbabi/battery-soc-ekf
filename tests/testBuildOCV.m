function tests = testBuildOCV
%TESTBUILDOCV Checks the approved degree-9 OCV fit.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.projectRoot = projectRoot;
addpath(fullfile(projectRoot, 'src'));
end

function testDegree9FitAndAnalyticDerivative(testCase)
[socTarget, voltageTarget] = approvedOCVTarget(testCase.TestData.projectRoot);
model = buildOCV(socTarget, voltageTarget);

denseSOC = linspace(0, 1, 10001)';
denseOCV = model.evaluate(denseSOC);
denseDerivative = model.evaluateDerivative(denseSOC);
targetFit = model.evaluate(socTarget);
targetError = targetFit - voltageTarget;

verifyEqual(testCase, model.degree, 9);
verifySize(testCase, model.coefficients, [1, 10]);
verifyEqual(testCase, model.derivativeCoefficients, ...
    polyder(model.coefficients));
verifySize(testCase, model.derivativeCoefficients, [1, 9]);
verifyTrue(testCase, all(isfinite(denseOCV)));
verifyTrue(testCase, all(isfinite(denseDerivative)));
verifyGreaterThan(testCase, min(denseDerivative), 0);
verifyLessThanOrEqual(testCase, model.fitQuality.rmse_V, 0.010);
verifyLessThanOrEqual(testCase, model.fitQuality.maxAbsError_V, 0.025);
verifyEqual(testCase, model.fitQuality.rmse_V, ...
    sqrt(mean(targetError .^ 2)), 'AbsTol', 1e-14);
verifyEqual(testCase, model.fitQuality.maxAbsError_V, ...
    max(abs(targetError)), 'AbsTol', 1e-14);
verifyEqual(testCase, model.fitQuality.endpointError_V, ...
    [targetError(1); targetError(end)], 'AbsTol', 1e-14);

fprintf('\nDegree-9 OCV fit evidence\n');
fprintf('  Target-grid RMSE: %.9f V\n', model.fitQuality.rmse_V);
fprintf('  Maximum absolute target error: %.9f V\n', ...
    model.fitQuality.maxAbsError_V);
fprintf('  SOC 0 endpoint error: %+.9f V\n', ...
    model.fitQuality.endpointError_V(1));
fprintf('  SOC 1 endpoint error: %+.9f V\n', ...
    model.fitQuality.endpointError_V(2));
fprintf('  Minimum dense-grid derivative: %.9f V per unit SOC\n', ...
    min(denseDerivative));
end

function [socTarget, voltageTarget] = approvedOCVTarget(projectRoot)
% Just calls into src/buildOCVTarget.m now (pulled out 2026-09-17 so this
% wasn't duplicated between the test and the real code). Same logic as before.
[socTarget, voltageTarget] = buildOCVTarget(string(projectRoot));
end
