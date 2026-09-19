function tests = testEcmSimulate
%TESTECMSIMULATE Right-endpoint 1-RC ECM contract tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
end

function testNonuniformRightEndpointRecurrenceAndVoltage(testCase)
t_s = [0; 0.7; 2.1; 3.0; 5.6; 7.0];
I_A = [0; 0; 2; 2; 0; 0];
SOC_ref = [0.80; 0.80; 0.772; 0.754; 0.754; 0.754];
params = syntheticParams();
evaluateOCV = @(soc) 3.20 + 0.50 * soc;

[V_model, V1] = ecmSimulate(t_s, I_A, SOC_ref, evaluateOCV, params);

expectedV1 = zeros(size(t_s));
for k = 2:numel(t_s)
    dt_s = t_s(k) - t_s(k - 1);
    alpha = exp(-dt_s / params.tau);
    expectedV1(k) = alpha * expectedV1(k - 1) + ...
        (1 - alpha) * params.R1 * I_A(k);
end
expectedVoltage = nan(size(t_s));
expectedVoltage(2:end) = evaluateOCV(SOC_ref(2:end)) - ...
    params.R0 * I_A(2:end) - expectedV1(2:end);

verifyEqual(testCase, V1, expectedV1, 'AbsTol', 1e-14);
verifyTrue(testCase, isnan(V_model(1)));
verifyEqual(testCase, V_model(2:end), expectedVoltage(2:end), 'AbsTol', 1e-14);
verifySize(testCase, V1, size(t_s));
verifySize(testCase, V_model, size(t_s));
verifyEqual(testCase, SOC_ref, [0.80; 0.80; 0.772; 0.754; 0.754; 0.754]);
end

function testRestToPulseUsesRightEndpointCurrent(testCase)
t_s = [0; 0.7; 2.1; 3.0];
I_A = [0; 0; 2; 2];
SOC_ref = [0.80; 0.80; 0.772; 0.754];
params = syntheticParams();
evaluateOCV = @(soc) 3.20 + 0.50 * soc;

[V_model, V1] = ecmSimulate(t_s, I_A, SOC_ref, evaluateOCV, params);

pulseStart = 3;
alphaAtPulse = exp(-(t_s(pulseStart) - t_s(pulseStart - 1)) / params.tau);
expectedV1AtPulse = (1 - alphaAtPulse) * params.R1 * I_A(pulseStart);
expectedVoltageAtPulse = evaluateOCV(SOC_ref(pulseStart)) - ...
    params.R0 * I_A(pulseStart) - expectedV1AtPulse;

verifyGreaterThan(testCase, V1(pulseStart), 0);
verifyEqual(testCase, V1(pulseStart), expectedV1AtPulse, 'AbsTol', 1e-14);
verifyEqual(testCase, V_model(pulseStart), expectedVoltageAtPulse, ...
    'AbsTol', 1e-14);
end

function testFirstSamplePolicyIsExplicitAndIndependent(testCase)
t_s = [0; 1.3];
I_A = [5; 0];
SOC_ref = [0.90; 0.80];
params = syntheticParams();

[V_model, V1] = ecmSimulate(t_s, I_A, SOC_ref, ...
    @ocvAfterFirstSampleOnly, params);

verifyEqual(testCase, V1(1), 0);
verifyTrue(testCase, isnan(V_model(1)));
verifyEqual(testCase, V_model(2), 3.60, 'AbsTol', 1e-14);
end

function testInputGuardErrorIdentifiers(testCase)
params = syntheticParams();
validTime_s = [0; 1; 2];
validCurrent_A = [0; 1; 0];
validSOC = [0.80; 0.79; 0.79];

verifyError(testCase, ...
    @() ecmSimulate(zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    @linearOCV, params), "ecmSimulate:EmptyInput");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A(1:2), validSOC, ...
    @linearOCV, params), "ecmSimulate:SizeMismatch");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, [0; NaN; 0], validSOC, ...
    @linearOCV, params), "ecmSimulate:NonFiniteVector");
verifyError(testCase, ...
    @() ecmSimulate([0; 0; 2], validCurrent_A, validSOC, ...
    @linearOCV, params), "ecmSimulate:InvalidTime");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, [0.80; -0.01; 0.79], ...
    @ocvMustNotRun, params), "ecmSimulate:InvalidReferenceSOC");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, [0.80; 1.01; 0.79], ...
    @ocvMustNotRun, params), "ecmSimulate:InvalidReferenceSOC");

missingTau = params;
missingTau = rmfield(missingTau, 'tau');
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, validSOC, ...
    @linearOCV, missingTau), "ecmSimulate:MissingParameter");

invalidParameter = params;
invalidParameter.R1 = 0;
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, validSOC, ...
    @linearOCV, invalidParameter), "ecmSimulate:InvalidParameter");

verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, validSOC, ...
    @wrongSizeOCV, params), "ecmSimulate:InvalidOCVOutput");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, validSOC, ...
    @nonfiniteOCV, params), "ecmSimulate:InvalidOCVOutput");
verifyError(testCase, ...
    @() ecmSimulate(validTime_s, validCurrent_A, validSOC, ...
    @nonvectorOCV, params), "ecmSimulate:InvalidOCVOutput");
end

function params = syntheticParams()
params = struct;
params.R0 = 0.010;
params.R1 = 0.040;
params.tau = 2.5;
end

function voltage_V = ocvAfterFirstSampleOnly(soc)
assert(all(soc < 0.90), ...
    'The OCV evaluator must not receive the first-sample SOC.');
voltage_V = 3.20 + 0.50 * soc;
end

function voltage_V = linearOCV(soc)
voltage_V = 3.20 + 0.50 * soc;
end

function voltage_V = ocvMustNotRun(~)
error('testEcmSimulate:OCVCalled', ...
    'The OCV evaluator must not run for invalid reference SOC.');
end

function voltage_V = wrongSizeOCV(~)
voltage_V = [3.4; 3.5; 3.6];
end

function voltage_V = nonfiniteOCV(soc)
voltage_V = nan(size(soc));
end

function voltage_V = nonvectorOCV(soc)
voltage_V = ones(numel(soc), 2);
end
