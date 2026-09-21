%EXPB_INITSOC Experiment B - wrong initial SOC: baseline SOC0 = 0.8 fed to
%both SOC_CC and SOC_EKF, while SOC_ref is built separately from the
%correct protocol SOC0 = 1. Current and voltage are otherwise clean (same
%conditions as Experiment A, so any difference from Experiment A's result
%is attributable only to the wrong initial condition).
%
%   Question this experiment answers: does the EKF recover from a wrong
%   initial SOC guess, and does Coulomb counting? See
%   tests/checkUS06CoulombCountInitialSOC.m for the original CC-only
%   evidence (constant -0.2 offset, exact to 5.55e-17) that this script
%   extends with the EKF side. See docs/DECISIONS.md, "Experiment B"
%   entry, for how the result below is read.

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

paramsHPPC = ecmParamsHPPC();

%% Reference SOC: correct protocol initial condition, clean current
SOC_ref = referenceSOC(data.t, data.I, 1.0, paramsHPPC.Q_As);

%% Wrong initial SOC fed to both baselines (0.8 instead of the true 1.0)
SOC0_wrong = 0.8;
SOC_CC = coulombCount(data.t, data.I, SOC0_wrong, paramsHPPC.Q_As);

%% EKF with the same wrong initial SOC. P0_SOC reflects that this initial
%  guess is known to be uncertain (20 percentage points), so the filter is
%  told up front not to trust it too much. That's what lets the voltage
%  measurements pull SOC back toward truth over time.
tuning = struct;
tuning.P0_SOC = (0.20)^2;  % large initial uncertainty: SOC0 is known-suspect here
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;
tuning.q_V1   = 1e-6;
tuning.R_V    = (3e-3)^2;

ekf = ekfEstimate(data.t, data.I, data.V, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsHPPC, SOC0_wrong, paramsHPPC.Q_As, tuning);

%% Metrics (whole-record, for comparability with Experiment A)
m_ekf = metrics(ekf.SOC, SOC_ref);
m_cc  = metrics(SOC_CC, SOC_ref);

err_ekf = ekf.SOC - SOC_ref;
err_cc  = SOC_CC - SOC_ref;
coverage3Sigma = 100 * mean(abs(err_ekf) <= ekf.sigma3_SOC);

% Convergence check: does the EKF's error shrink toward Experiment A's
% steady-state error level, or stay pinned near the initial -0.2 offset?
initialErr = err_ekf(1);
lateWindow = round(0.5 * numel(data.t)):numel(data.t);
lateMeanAbsErr = mean(abs(err_ekf(lateWindow)));

fprintf('\nExperiment B - wrong initial SOC (SOC0 = %.2f fed to CC and EKF, true SOC0 = 1.0)\n', SOC0_wrong);
fprintf('  EKF  RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_ekf.RMSE, m_ekf.MAE, m_ekf.bias, m_ekf.maxAbsError);
fprintf('  CC   RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_cc.RMSE, m_cc.MAE, m_cc.bias, m_cc.maxAbsError);
fprintf('  EKF +/-3 sigma coverage: %.4f %%\n', coverage3Sigma);
fprintf('  EKF initial error: %+.6f, EKF mean |error| over the back half of the record: %.6f\n', ...
    initialErr, lateMeanAbsErr);
fprintf('  CC error is constant by construction: initial = %+.6f, final = %+.6f\n', ...
    err_cc(1), err_cc(end));

fprintf(['\nInterpretation (see DECISIONS.md): CC''s initial and final ' ...
    'error match to numerical precision, confirming it has no way to ' ...
    'correct a wrong initial condition. Once it starts -0.2 off, it stays ' ...
    '-0.2 off for the whole record, and that constant offset is basically ' ...
    'what CC''s RMSE here is measuring. The EKF is different: the Kalman ' ...
    'gain keeps pulling the estimate toward whatever the voltage implies ' ...
    'SOC should be, so the error shrinks from that initial -0.2 down ' ...
    'toward Experiment A''s steady-state level as the record goes on. ' ...
    'That is the real advantage an estimator has over plain integration. ' ...
    'It does not reach zero, though; it settles toward Experiment A''s ' ...
    'bias level, because the same ~72 mV HPPC-ECM voltage bias from that ' ...
    'experiment is still baked in here and puts a floor on how precisely ' ...
    'voltage can pin down SOC.\n']);

%% Save results and figure
if ~exist(fullfile(projectRoot, 'results'), 'dir')
    mkdir(fullfile(projectRoot, 'results'));
end
save(fullfile(projectRoot, 'results', 'expB.mat'), ...
    'SOC_ref', 'SOC_CC', 'ekf', 'm_ekf', 'm_cc', 'coverage3Sigma', 'SOC0_wrong');

figure('Position', [100, 100, 900, 600], 'Color', 'w');
try, set(gcf, 'Theme', 'light'); end  %#ok<TRYNC>  % keep export light even if the MATLAB desktop is in dark mode
subplot(2, 1, 1);
plot(data.t, SOC_ref, 'k-', 'LineWidth', 1.5); hold on;
plot(data.t, ekf.SOC, 'b-');
plot(data.t, SOC_CC, 'r--');
xlabel('Time (s)'); ylabel('SOC (-)');
legend('Reference', 'EKF (wrong SOC0)', 'Coulomb counting (wrong SOC0)', 'Location', 'best');
grid on; title(sprintf('Experiment B - wrong initial SOC (SOC0 = %.2f, true = 1.0)', SOC0_wrong));

subplot(2, 1, 2);
plot(data.t, err_ekf, 'b-'); hold on;
plot(data.t, err_cc, 'r--');
yline(0, 'k:');
% CC's error sits at exactly -0.2 for the whole record, so without this the
% line is drawn along the axis frame and is effectively invisible.
ylim([-0.24, 0.05]);
xlabel('Time (s)'); ylabel('SOC error (-)');
legend('EKF error', 'CC error (constant offset)', 'Location', 'best');
grid on; title('EKF recovers from the wrong initial condition; CC does not');

if ~exist(fullfile(projectRoot, 'figures'), 'dir')
    mkdir(fullfile(projectRoot, 'figures'));
end
exportgraphics(gcf, fullfile(projectRoot, 'figures', 'expB_initSOC.png'), 'Resolution', 300, 'BackgroundColor', 'white');
