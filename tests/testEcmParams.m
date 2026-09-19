function tests = testEcmParams
%TESTECMPARAMS Checks the approved provisional ECM constants.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
testCase.TestData.expectedTau = 6.343342839558;
testCase.TestData.expectedQAs = 10790.6154948;
end

function testExactApprovedValuesAndDerivedRelations(testCase)
params = ecmParams();

verifyEqual(testCase, params.R0, 0.004500052);
verifyEqual(testCase, params.R1, 0.028981830);
verifyEqual(testCase, params.C1, 218.8730953);
verifyEqual(testCase, params.Q_Ah, 2.997393193);
verifyEqual(testCase, params.tau, testCase.TestData.expectedTau, ...
    'AbsTol', 5e-13);
verifyEqual(testCase, params.Q_As, testCase.TestData.expectedQAs, ...
    'AbsTol', 5e-10);

numericValues = [params.R0, params.R1, params.C1, params.tau, ...
    params.Q_Ah, params.Q_As];
verifyTrue(testCase, all(isfinite(numericValues)));
verifyTrue(testCase, all(numericValues > 0));
verifyEqual(testCase, params.tau, params.R1 * params.C1, ...
    'AbsTol', eps(params.tau));
verifyEqual(testCase, params.Q_As, params.Q_Ah * 3600, ...
    'AbsTol', eps(params.Q_As));
end
