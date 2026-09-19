%EXPA_NOMINAL Experiment A - nominal case: correct init SOC, HPPC params,
%no added corruption. Compares SOC_EKF and SOC_CC against SOC_ref on the
%approved US06 record.
%
%   This is the first EKF experiment in the project and the baseline every
%   later experiment (B-E) gets compared against. See docs/DECISIONS.md,
%   "Starting the EKF on an imperfect ECM" and "Experiment A nominal
%   result" for why the nominal parameter set is ecmParamsHPPC() instead
%   of ecmParams(), and for how to read the residual bias printed below
%   (it isn't tuned away here, on purpose).

clear; clc;
rng(42);

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

%% Load approved data and build the OCV model
us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[data, ~] = loadDataset(string(us06File));

[socTarget, voltageTarget] = buildOCVTarget(string(projectRoot));
ocvModel = buildOCV(socTarget, voltageTarget);

%% Reference SOC (clean current, protocol-supported SOC_ref(1) = 1)
paramsHPPC = ecmParamsHPPC();
SOC_ref = referenceSOC(data.t, data.I, 1.0, paramsHPPC.Q_As);

%% Coulomb-counting baseline (same clean current as Experiment A; no
%  corruption yet, that starts at Experiment B/C/E)
SOC_CC = coulombCount(data.t, data.I, 1.0, paramsHPPC.Q_As);

%% EKF estimate, nominal tuning
tuning = struct;
tuning.P0_SOC = (0.01)^2;   % 1% initial-SOC uncertainty (SOC0 is correct here)
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;      % SOC dynamics trusted almost exactly (clean current)
tuning.q_V1   = 1e-6;       % more model uncertainty in the RC branch
tuning.R_V    = (3e-3)^2;   % measurement-noise variance, 3 mV std

ekf = ekfEstimate(data.t, data.I, data.V, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsHPPC, 1.0, paramsHPPC.Q_As, tuning);

%% Metrics
m_ekf = metrics(ekf.SOC, SOC_ref);
m_cc  = metrics(SOC_CC, SOC_ref);

err_ekf = ekf.SOC - SOC_ref;
coverage3Sigma = 100 * mean(abs(err_ekf) <= ekf.sigma3_SOC);

fprintf('\nExperiment A - nominal (HPPC parameters, correct init, no added corruption)\n');
fprintf('  EKF  RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_ekf.RMSE, m_ekf.MAE, m_ekf.bias, m_ekf.maxAbsError);
fprintf('  CC   RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_cc.RMSE, m_cc.MAE, m_cc.bias, m_cc.maxAbsError);
fprintf('  EKF +/-3 sigma coverage: %.4f %%\n', coverage3Sigma);

fprintf(['\nInterpretation (see DECISIONS.md): SOC_CC''s zero error here ' ...
    'isn''t a real result. It''s built from the same clean current as ' ...
    'SOC_ref, so the two integrate identically by construction, and CC ' ...
    'only starts telling us something once its input is corrupted, ' ...
    'starting in Experiment B. The EKF''s bias comes from somewhere real: ' ...
    'the filter is absorbing the HPPC-identified ECM''s own residual ' ...
    'voltage bias (+0.0717 V, tests/checkUS06HPPC1RCVoltage.m) through ' ...
    'dOCV/dSOC, so this is expected given that the ECM fit has that bias ' ...
    'to begin with. The low +/-3 sigma coverage shows the filter is ' ...
    'overconfident under this process-model mismatch, and Q was left ' ...
    'alone here, not inflated to paper over it.\n']);

%% Save results and figure
if ~exist(fullfile(projectRoot, 'results'), 'dir')
    mkdir(fullfile(projectRoot, 'results'));
end
save(fullfile(projectRoot, 'results', 'expA.mat'), ...
    'SOC_ref', 'SOC_CC', 'ekf', 'm_ekf', 'm_cc', 'coverage3Sigma');

figure('Position', [100, 100, 900, 600]);
subplot(2, 1, 1);
plot(data.t, SOC_ref, 'k-', 'LineWidth', 1.5); hold on;
plot(data.t, ekf.SOC, 'b-');
plot(data.t, SOC_CC, 'r--');
xlabel('Time (s)'); ylabel('SOC (-)');
legend('Reference', 'EKF', 'Coulomb counting', 'Location', 'best');
grid on; title('Experiment A - nominal case');

subplot(2, 1, 2);
fill([data.t; flipud(data.t)], [ekf.sigma3_SOC; flipud(-ekf.sigma3_SOC)], ...
    [0.85 0.85 0.95], 'EdgeColor', 'none'); hold on;
plot(data.t, err_ekf, 'b-');
xlabel('Time (s)'); ylabel('SOC error (-)');
legend('+/-3\sigma band', 'EKF error', 'Location', 'best');
grid on; title(sprintf('EKF error vs +/-3\\sigma band (coverage = %.2f%%)', coverage3Sigma));

if ~exist(fullfile(projectRoot, 'figures'), 'dir')
    mkdir(fullfile(projectRoot, 'figures'));
end
exportgraphics(gcf, fullfile(projectRoot, 'figures', 'expA_nominal.png'), 'Resolution', 300);
