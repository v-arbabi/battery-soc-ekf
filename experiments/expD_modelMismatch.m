%EXPD_MODELMISMATCH Experiment D - process model swapped for a worse fit:
%the EKF is given the literature-derived 1-RC parameters (src/ecmParams.m,
%SOC = 0.5 constant set) instead of the HPPC-identified nominal set used
%in Experiment A. Current, voltage, and initial SOC are otherwise clean
%(identical inputs to Experiment A), so any change in EKF performance is
%attributable only to the process-model parameters.
%
%   Question this experiment answers: how much worse does the EKF get when
%   given the weaker-fitting parameter set? See docs/DECISIONS.md,
%   "EKF entry authorization" for why the literature set was earmarked for
%   this experiment, and the "Experiment D" entry for the result.

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

paramsHPPC = ecmParamsHPPC();          % used only for the shared reference capacity
paramsLiterature = ecmParams();        % the weaker-fitting model under test, SOC=0.5 constants

assert(abs(paramsHPPC.Q_As - paramsLiterature.Q_As) < 1e-6, ...
    'expD_modelMismatch:CapacityMismatch', ...
    'Both parameter sets must share the same measured capacity so this experiment isolates R0/R1/tau only.');

%% Reference SOC and Coulomb-counting baseline (clean current, correct SOC0)
SOC_ref = referenceSOC(data.t, data.I, 1.0, paramsHPPC.Q_As);
SOC_CC  = coulombCount(data.t, data.I, 1.0, paramsHPPC.Q_As);

%% EKF, given the mismatched (literature) process-model parameters
tuning = struct;
tuning.P0_SOC = (0.01)^2;
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;
tuning.q_V1   = 1e-6;
tuning.R_V    = (3e-3)^2;

ekf = ekfEstimate(data.t, data.I, data.V, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsLiterature, 1.0, paramsHPPC.Q_As, tuning);

%% Metrics, compared directly against Experiment A's numbers
m_ekf = metrics(ekf.SOC, SOC_ref);
m_cc  = metrics(SOC_CC, SOC_ref);

err_ekf = ekf.SOC - SOC_ref;
coverage3Sigma = 100 * mean(abs(err_ekf) <= ekf.sigma3_SOC);

fprintf('\nExperiment D - deliberately mismatched process model (literature R0/R1/tau vs. HPPC in Exp A)\n');
fprintf('  EKF  RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_ekf.RMSE, m_ekf.MAE, m_ekf.bias, m_ekf.maxAbsError);
fprintf('  CC   RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_cc.RMSE, m_cc.MAE, m_cc.bias, m_cc.maxAbsError);
fprintf('  EKF +/-3 sigma coverage: %.4f %%\n', coverage3Sigma);

fprintf(['\nInterpretation (see DECISIONS.md): the literature parameter ' ...
    'set fits the US06 voltage worse than the HPPC set from Experiment A ' ...
    'does: +0.089 V bias and 0.114 V RMSE, versus +0.072 V and 0.084 V ' ...
    'for HPPC (tests/checkUS06ECMVoltage.m vs. ' ...
    'tests/checkUS06HPPC1RCVoltage.m). ' ...
    'Everything else here (current, voltage, initial SOC) matches ' ...
    'Experiment A exactly, so this experiment isolates whether that ' ...
    'worse voltage fit shows up as worse EKF SOC accuracy. CC is not ' ...
    'affected by any of this: same clean current, same capacity as ' ...
    'Experiment A, so it is included only as the fixed reference point, ' ...
    'not as part of the comparison.\n']);

%% Save results and figure
if ~exist(fullfile(projectRoot, 'results'), 'dir')
    mkdir(fullfile(projectRoot, 'results'));
end
save(fullfile(projectRoot, 'results', 'expD.mat'), ...
    'SOC_ref', 'SOC_CC', 'ekf', 'm_ekf', 'm_cc', 'coverage3Sigma');

figure('Position', [100, 100, 900, 600]);
subplot(2, 1, 1);
plot(data.t, SOC_ref, 'k-', 'LineWidth', 1.5); hold on;
plot(data.t, ekf.SOC, 'b-');
plot(data.t, SOC_CC, 'r--');
xlabel('Time (s)'); ylabel('SOC (-)');
legend('Reference', 'EKF (literature params)', 'Coulomb counting', 'Location', 'best');
grid on; title('Experiment D - deliberately mismatched process model');

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
exportgraphics(gcf, fullfile(projectRoot, 'figures', 'expD_modelMismatch.png'), 'Resolution', 300);
