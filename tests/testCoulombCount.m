function tests = testCoulombCount
%TESTCOULOMBCOUNT Tests for the current-only SOC_CC integration.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
end

function testNonuniformRightEndpointIntegrationAndInputPreservation(testCase)
t_s = [0; 1.5; 2.0; 5.0];
I_A = [0; 2; -1; 3];
SOC0 = 0.80;
Q_As = 100;
originalTime_s = t_s;
originalCurrent_A = I_A;

SOC_CC = coulombCount(t_s, I_A, SOC0, Q_As);

expectedSOC = [0.80; 0.77; 0.775; 0.685];
verifyEqual(testCase, SOC_CC, expectedSOC, 'AbsTol', 1e-14);
verifySize(testCase, SOC_CC, size(t_s));
verifyEqual(testCase, SOC_CC(1), SOC0);
verifyLessThan(testCase, SOC_CC(2), SOC_CC(1));
verifyGreaterThan(testCase, SOC_CC(3), SOC_CC(2));
verifyEqual(testCase, t_s, originalTime_s);
verifyEqual(testCase, I_A, originalCurrent_A);
end

function testOutputIsNotClamped(testCase)
SOC_CC = coulombCount([0; 2], [0; 1], 0.5, 1);

verifyEqual(testCase, SOC_CC, [0.5; -1.5], 'AbsTol', 1e-14);
verifyLessThan(testCase, SOC_CC(end), 0);
end

function testInputValidationErrorIdentifiers(testCase)
verifyError(testCase, ...
    @() coulombCount(zeros(0, 1), zeros(0, 1), 0.5, 100), ...
    "coulombCount:EmptyInput");
verifyError(testCase, ...
    @() coulombCount([0; 1], 0, 0.5, 100), ...
    "coulombCount:SizeMismatch");
verifyError(testCase, ...
    @() coulombCount([0; 1], [0; NaN], 0.5, 100), ...
    "coulombCount:NonFiniteVector");
verifyError(testCase, ...
    @() coulombCount([0; 0], [0; 1], 0.5, 100), ...
    "coulombCount:InvalidTime");
verifyError(testCase, ...
    @() coulombCount([0; 1], [0; 1], 1.1, 100), ...
    "coulombCount:InvalidInitialSOC");
verifyError(testCase, ...
    @() coulombCount([0; 1], [0; 1], 0.5, 0), ...
    "coulombCount:InvalidCapacity");
end
