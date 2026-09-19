%EXPC_SENSOR Experiment C - corrupted sensors: a constant +0.02 A current
%bias (same magnitude as tests/checkUS06CoulombCountCurrentBias.m) plus
%zero-mean Gaussian voltage measurement noise (3 mV std, matching the
%EKF's own R_V assumption) are fed to both SOC_CC and SOC_EKF. SOC_ref is
%still built from the clean current, so it is unaffected by either
%corruption and stays a fair comparison target.
%
%   Question this experiment answers: with a biased current sensor and a
%   noisy voltage sensor, does fusing the (biased) current with the
%   (noisy but unbiased-in-expectation) voltage measurement do better than
%   Coulomb counting on the biased current alone? See docs/DECISIONS.md,
%   "Experiment C" entry for the result and why it isn't what I first
%   expected.

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

%% Reference SOC: clean current, correct initial SOC (unaffected by the
%  sensor corruption introduced below)
SOC_ref = referenceSOC(data.t, data.I, 1.0, paramsHPPC.Q_As);

%% Corrupted sensor signals
currentBias_A = 0.02;      % same fixed assumption as Experiment-B-era CC evidence
voltageNoiseStd_V = 3e-3;  % matches tuning.R_V below; explicit assumption, not measured

I_corrupted = data.I + currentBias_A;
V_corrupted = data.V + voltageNoiseStd_V * randn(size(data.V));

%% Coulomb-counting baseline sees only the biased current
SOC_CC = coulombCount(data.t, I_corrupted, 1.0, paramsHPPC.Q_As);

%% EKF sees both the biased current and the noisy voltage
tuning = struct;
tuning.P0_SOC = (0.01)^2;
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;
tuning.q_V1   = 1e-6;
tuning.R_V    = (voltageNoiseStd_V)^2;  % filter is told the true noise level

ekf = ekfEstimate(data.t, I_corrupted, V_corrupted, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsHPPC, 1.0, paramsHPPC.Q_As, tuning);

%% Metrics
m_ekf = metrics(ekf.SOC, SOC_ref);
m_cc  = metrics(SOC_CC, SOC_ref);

err_ekf = ekf.SOC - SOC_ref;
err_cc  = SOC_CC - SOC_ref;
coverage3Sigma = 100 * mean(abs(err_ekf) <= ekf.sigma3_SOC);

fprintf('\nExperiment C - corrupted sensors (current bias +%.3f A, voltage noise %.1f mV std)\n', ...
    currentBias_A, voltageNoiseStd_V * 1000);
fprintf('  EKF  RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_ekf.RMSE, m_ekf.MAE, m_ekf.bias, m_ekf.maxAbsError);
fprintf('  CC   RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_cc.RMSE, m_cc.MAE, m_cc.bias, m_cc.maxAbsError);
fprintf('  EKF +/-3 sigma coverage: %.4f %%\n', coverage3Sigma);

fprintf(['\nInterpretation (see DECISIONS.md): CC comes out ahead in this ' ...
    'one. Its error stays under 1%% because a +0.02 A current bias just ' ...
    'isn''t very large, while the EKF''s error barely moves from ' ...
    'Experiment A (still around 7%%), which tells you what is actually ' ...
    'driving it: the same ~72 mV HPPC-ECM voltage-fit bias, not the ' ...
    'sensor corruption added in this experiment. I expected fusing in the ' ...
    'voltage measurement to win out over plain integration here, but that ' ...
    'is not what happened, because the corruption injected is small next ' ...
    'to the model bias the filter was already carrying before this ' ...
    'experiment started. Keeping this as a CC win rather than dressing it ' ...
    'up as an EKF success is the point: it is a genuine case where sensor ' ...
    'fusion does not help, and it belongs in the report as-is.\n']);

%% Save results and figure
if ~exist(fullfile(projectRoot, 'results'), 'dir')
    mkdir(fullfile(projectRoot, 'results'));
end
save(fullfile(projectRoot, 'results', 'expC.mat'), ...
    'SOC_ref', 'SOC_CC', 'ekf', 'm_ekf', 'm_cc', 'coverage3Sigma', ...
    'currentBias_A', 'voltageNoiseStd_V');

figure('Position', [100, 100, 900, 600], 'Color', 'w');
try, set(gcf, 'Theme', 'light'); end  %#ok<TRYNC>  % keep export light even if the MATLAB desktop is in dark mode
subplot(2, 1, 1);
plot(data.t, SOC_ref, 'k-', 'LineWidth', 1.5); hold on;
plot(data.t, ekf.SOC, 'b-');
plot(data.t, SOC_CC, 'r--');
xlabel('Time (s)'); ylabel('SOC (-)');
legend('Reference', 'EKF (corrupted sensors)', 'Coulomb counting (biased current)', 'Location', 'best');
grid on; title('Experiment C - corrupted current + voltage sensors');

subplot(2, 1, 2);
plot(data.t, err_ekf, 'b-'); hold on;
plot(data.t, err_cc, 'r--');
yline(0, 'k:');
xlabel('Time (s)'); ylabel('SOC error (-)');
legend('EKF error', 'CC error', 'Location', 'best');
grid on; title('CC error stays small here; EKF error is dominated by pre-existing model bias, not by this sensor fault');

if ~exist(fullfile(projectRoot, 'figures'), 'dir')
    mkdir(fullfile(projectRoot, 'figures'));
end
exportgraphics(gcf, fullfile(projectRoot, 'figures', 'expC_sensor.png'), 'Resolution', 300, 'BackgroundColor', 'white');
