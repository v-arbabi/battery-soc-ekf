function tests = testReferenceSOC
%TESTREFERENCESOC Actual-time reference SOC integration tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.projectRoot = projectRoot;
addpath(fullfile(projectRoot, 'src'));
end

function testNonuniformSyntheticRightEndpointIntegration(testCase)
t = [0; 0.5; 2.0; 5.5];
I_clean = [0; 2; -1; 4];
SOC0 = 0.8;
Q_As = 100;
expectedSOC = [0.8; 0.79; 0.805; 0.665];

SOC_ref = referenceSOC(t, I_clean, SOC0, Q_As);

verifyEqual(testCase, SOC_ref, expectedSOC, 'AbsTol', 2e-16);
verifyEqual(testCase, SOC_ref(1), SOC0);

% referenceSOC shouldn't clamp; check it still returns an out-of-range
% value when the inputs actually call for one.
unclampedSOC = referenceSOC([0; 2], [0; 1], 0.5, 1);
verifyLessThan(testCase, unclampedSOC(end), 0);
end

function testFullRetainedUS06Trajectory(testCase)
us06File = fullfile(testCase.TestData.projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[data, report] = loadDataset(string(us06File));
params = ecmParams();

SOC_ref = referenceSOC(data.t, data.I, 1, params.Q_As);

verifyEqual(testCase, numel(SOC_ref), 48060);
verifyEqual(testCase, numel(SOC_ref), numel(data.t));
verifyEqual(testCase, report.removedDuplicateSourceRows, 48061);
verifyEqual(testCase, SOC_ref(1), 1);
verifyTrue(testCase, all(isfinite(SOC_ref)));
verifyEqual(testCase, SOC_ref(end), 0.137215622136, ...
    'AbsTol', 5e-13);

expectedFromActualIntervals = 1 - ...
    sum(data.I(2:end) .* diff(data.t)) / params.Q_As;
verifyEqual(testCase, SOC_ref(end), expectedFromActualIntervals, ...
    'AbsTol', 1e-14);

fprintf('\nUS06 reference SOC evidence\n');
fprintf('  Retained samples: %d\n', numel(SOC_ref));
fprintf('  Initial SOC_ref: %.12f\n', SOC_ref(1));
fprintf('  Final SOC_ref: %.12f\n', SOC_ref(end));
end

function testInputValidationErrorIdentifiers(testCase)
verifyError(testCase, ...
    @() referenceSOC(zeros(0, 1), zeros(0, 1), 0.5, 100), ...
    "referenceSOC:EmptyInput");
verifyError(testCase, ...
    @() referenceSOC([0; 1], 0, 0.5, 100), ...
    "referenceSOC:SizeMismatch");
verifyError(testCase, ...
    @() referenceSOC([0; 1], [0; NaN], 0.5, 100), ...
    "referenceSOC:NonFiniteVector");
verifyError(testCase, ...
    @() referenceSOC([0; 0], [0; 1], 0.5, 100), ...
    "referenceSOC:InvalidTime");
verifyError(testCase, ...
    @() referenceSOC([0; 1], [0; 1], 1.1, 100), ...
    "referenceSOC:InvalidInitialSOC");
verifyError(testCase, ...
    @() referenceSOC([0; 1], [0; 1], 0.5, 0), ...
    "referenceSOC:InvalidCapacity");
end
