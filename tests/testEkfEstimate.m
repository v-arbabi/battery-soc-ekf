function tests = testEkfEstimate
%TESTEKFESTIMATE Tests for the two-state EKF SOC estimator.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
end

function testFirstSampleInitialisation(testCase)
t_s = [0; 1; 2];
I_A = [0; 1; 1];
V_meas = [3.7; 3.65; 3.6];
evaluateOCV = @(soc) 3.0 + 1.0 * soc;
evaluateDOCV = @(soc) 1.0 + 0 * soc;
params = syntheticParams();
tuning = syntheticTuning();

result = ekfEstimate(t_s, I_A, V_meas, evaluateOCV, evaluateDOCV, ...
    params, 0.6, 3600, tuning);

verifyEqual(testCase, result.SOC(1), 0.6);
verifyEqual(testCase, result.V1(1), 0);
verifyEqual(testCase, result.P_SOC(1), tuning.P0_SOC);
verifyTrue(testCase, isnan(result.innovation(1)));
verifySize(testCase, result.SOC, size(t_s));
verifySize(testCase, result.sigma3_SOC, size(t_s));
verifyEqual(testCase, result.sigma3_SOC, 3 * sqrt(result.P_SOC), 'AbsTol', 1e-14);
end

function testRecurrenceMatchesIndependentReferenceComputation(testCase)
t_s = [0; 0.5; 1.3; 2.0; 4.5; 5.0];
I_A = [0; 1.5; 1.5; -0.5; 0.8; 0.8];
params = syntheticParams();
evaluateOCV = @(soc) 3.0 + 1.2 * soc;
evaluateDOCV = @(soc) 1.2 + 0 * soc;
tuning = syntheticTuning();
SOC0 = 0.75;
Q_As = 3600;

result = ekfEstimate(t_s, I_A, [3.72; 3.70; 3.65; 3.66; 3.55; 3.54], ...
    evaluateOCV, evaluateDOCV, params, SOC0, Q_As, tuning);

% Reference computation done by hand, straight from the docstring
% equations, instead of reaching into ekfEstimate's internals.
V_meas = [3.72; 3.70; 3.65; 3.66; 3.55; 3.54];
x = [SOC0; 0];
P = diag([tuning.P0_SOC, tuning.P0_V1]);
Qk = diag([tuning.q_SOC, tuning.q_V1]);
Rk = tuning.R_V;
N = numel(t_s);
expectedSOC = zeros(N, 1);
expectedSOC(1) = SOC0;
for k = 2:N
    dt_k = t_s(k) - t_s(k - 1);
    alpha = exp(-dt_k / params.tau);
    x_pred = [x(1) - (dt_k / Q_As) * I_A(k); ...
        alpha * x(2) + (1 - alpha) * params.R1 * I_A(k)];
    F = [1, 0; 0, alpha];
    P_pred = F * P * F' + Qk;
    y_pred = evaluateOCV(x_pred(1)) - x_pred(2) - params.R0 * I_A(k);
    resid = V_meas(k) - y_pred;
    H = [evaluateDOCV(x_pred(1)), -1];
    S = H * P_pred * H' + Rk;
    K = (P_pred * H') / S;
    x_upd = x_pred + K * resid;
    P = (eye(2) - K * H) * P_pred;
    x = x_upd;
    x(1) = min(max(x(1), 0), 1);
    expectedSOC(k) = x(1);
end

verifyEqual(testCase, result.SOC, expectedSOC, 'AbsTol', 1e-12);
end

function testSOCStaysWithinPhysicalBounds(testCase)
% Big forced discharge current: SOC should clamp at 0 instead of going negative.
t_s = (0:1:50)';
I_A = 50 * ones(size(t_s));
V_meas = 3.0 * ones(size(t_s));
evaluateOCV = @(soc) 3.0 + 1.0 * soc;
evaluateDOCV = @(soc) 1.0 + 0 * soc;
params = syntheticParams();
tuning = syntheticTuning();

result = ekfEstimate(t_s, I_A, V_meas, evaluateOCV, evaluateDOCV, ...
    params, 0.5, 3600, tuning);

verifyTrue(testCase, all(result.SOC >= 0 & result.SOC <= 1));
end

function testInputGuards(testCase)
params = syntheticParams();
tuning = syntheticTuning();
evaluateOCV = @(soc) 3.0 + soc;
evaluateDOCV = @(soc) 1 + 0 * soc;

verifyError(testCase, @() ekfEstimate([], [], [], evaluateOCV, ...
    evaluateDOCV, params, 0.5, 3600, tuning), "ekfEstimate:EmptyInput");
verifyError(testCase, @() ekfEstimate([0; 1], [0], [3.7; 3.6], ...
    evaluateOCV, evaluateDOCV, params, 0.5, 3600, tuning), ...
    "ekfEstimate:SizeMismatch");
verifyError(testCase, @() ekfEstimate([0; 1; 0.5], [0; 1; 1], ...
    [3.7; 3.6; 3.6], evaluateOCV, evaluateDOCV, params, 0.5, 3600, tuning), ...
    "ekfEstimate:InvalidTime");
verifyError(testCase, @() ekfEstimate([0; 1], [0; 1], [3.7; 3.6], ...
    evaluateOCV, evaluateDOCV, params, 1.5, 3600, tuning), ...
    "ekfEstimate:InvalidInitialSOC");
verifyError(testCase, @() ekfEstimate([0; 1], [0; 1], [3.7; 3.6], ...
    evaluateOCV, evaluateDOCV, params, 0.5, -1, tuning), ...
    "ekfEstimate:InvalidCapacity");
incompleteTuning = rmfield(tuning, 'R_V');
verifyError(testCase, @() ekfEstimate([0; 1], [0; 1], [3.7; 3.6], ...
    evaluateOCV, evaluateDOCV, params, 0.5, 3600, incompleteTuning), ...
    "ekfEstimate:MissingTuning");
end

function params = syntheticParams()
params = struct;
params.R0 = 0.02;
params.R1 = 0.015;
params.tau = 8;
end

function tuning = syntheticTuning()
tuning = struct;
tuning.P0_SOC = (0.05)^2;
tuning.P0_V1 = (0.01)^2;
tuning.q_SOC = 1e-9;
tuning.q_V1 = 1e-6;
tuning.R_V = (5e-3)^2;
end
