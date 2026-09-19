function tests = testMetrics
%TESTMETRICS Checks for the generic aligned-vector error stats.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
end

function testKnownStatisticsAndInputPreservation(testCase)
estimate = [0.9; 0.2; 1.0; -0.2];
reference = [0.7; 0.3; 0.5; 0.1];
estimateBefore = estimate;
referenceBefore = reference;

result = metrics(estimate, reference);
expectedError = [0.2; -0.1; 0.5; -0.3];

verifyEqual(testCase, result.error, expectedError, 'AbsTol', 1e-14);
verifyEqual(testCase, result.RMSE, sqrt(mean(expectedError .^ 2)), 'AbsTol', 1e-14);
verifyEqual(testCase, result.MAE, 0.275, 'AbsTol', 1e-14);
verifyEqual(testCase, result.bias, 0.075, 'AbsTol', 1e-14);
verifyEqual(testCase, result.maxAbsError, 0.5, 'AbsTol', 1e-14);
verifyEqual(testCase, result.validSampleCount, 4);
verifyEqual(testCase, estimate, estimateBefore);
verifyEqual(testCase, reference, referenceBefore);
end

function testRowAndColumnVectorsRemainAlignedByElementOrder(testCase)
estimate = [0.4, 0.1, 0.9];
reference = [0.3; 0.3; 0.2];

result = metrics(estimate, reference);

verifyEqual(testCase, result.error, [0.1; -0.2; 0.7], 'AbsTol', 1e-14);
verifySize(testCase, result.error, [3, 1]);
verifyEqual(testCase, result.validSampleCount, 3);
end

function testInputErrors(testCase)
validEstimate = [0.2; 0.4];
validReference = [0.1; 0.3];

verifyError(testCase, @() metrics(zeros(0, 1), zeros(0, 1)), ...
    "metrics:EmptyInput");
verifyError(testCase, @() metrics(validEstimate, validReference(1)), ...
    "metrics:SizeMismatch");
verifyError(testCase, @() metrics([0.2; NaN], validReference), ...
    "metrics:NonFiniteInput");
verifyError(testCase, @() metrics([0.2 + 0.1i; 0.4], validReference), ...
    "metrics:NonFiniteInput");
verifyError(testCase, @() metrics([0.2; Inf], validReference), ...
    "metrics:NonFiniteInput");
verifyError(testCase, @() metrics([0.2, 0.4; 0.3, 0.5], ...
    [0.1; 0.3; 0.2; 0.4]), "metrics:NonVectorInput");
end
