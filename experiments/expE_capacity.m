%EXPE_CAPACITY Experiment E - capacity mismatch: both SOC_CC and the EKF
%are given the documented nominal capacity (2.9 Ah) instead of the
%approved measured capacity (2.997393193 Ah) used to build SOC_ref. This
%extends tests/checkUS06CoulombCountCapacityMismatch.m (CC-only evidence)
%with the EKF side. Current, voltage, and initial SOC are otherwise clean.
%
%   Question this experiment answers: a wrong capacity biases the process
%   model's SOC integration step itself (dt/Q_As), not just a single
%   parameter. Does the EKF's voltage correction help here the way it did
%   in Experiment C, or does a wrong Q_As resist correction differently?
%   See docs/DECISIONS.md, "Experiment E" entry, for the result.

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

Q_As_measured = paramsHPPC.Q_As;         % 2.997393193 Ah, used for SOC_ref only
Q_Ah_nominal  = 2.9;                     % documented nominal, used for CC and EKF only
Q_As_nominal  = Q_Ah_nominal * 3600;

%% Reference SOC: correct measured capacity, clean current
SOC_ref = referenceSOC(data.t, data.I, 1.0, Q_As_measured);

%% Coulomb-counting baseline: same clean current, wrong (nominal) capacity
SOC_CC = coulombCount(data.t, data.I, 1.0, Q_As_nominal);

%% EKF: same clean current and voltage, wrong (nominal) capacity in its
%  own process model (Q_As enters the SOC prediction step directly)
tuning = struct;
tuning.P0_SOC = (0.01)^2;
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;
tuning.q_V1   = 1e-6;
tuning.R_V    = (3e-3)^2;

ekf = ekfEstimate(data.t, data.I, data.V, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsHPPC, 1.0, Q_As_nominal, tuning);

%% Metrics
m_ekf = metrics(ekf.SOC, SOC_ref);
m_cc  = metrics(SOC_CC, SOC_ref);

err_ekf = ekf.SOC - SOC_ref;
err_cc  = SOC_CC - SOC_ref;
coverage3Sigma = 100 * mean(abs(err_ekf) <= ekf.sigma3_SOC);

fprintf('\nExperiment E - capacity mismatch (CC and EKF given nominal %.2f Ah vs. measured %.9f Ah)\n', ...
    Q_Ah_nominal, paramsHPPC.Q_Ah);
fprintf('  EKF  RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_ekf.RMSE, m_ekf.MAE, m_ekf.bias, m_ekf.maxAbsError);
fprintf('  CC   RMSE = %.6f  MAE = %.6f  bias = %+.6f  maxAbsError = %.6f\n', ...
    m_cc.RMSE, m_cc.MAE, m_cc.bias, m_cc.maxAbsError);
fprintf('  EKF +/-3 sigma coverage: %.4f %%\n', coverage3Sigma);

fprintf(['\nInterpretation (see DECISIONS.md): a wrong Q_As scales every ' ...
    'prediction step''s SOC increment, (dt/Q_As)*I, so CC''s error grows ' ...
    'with how much current has flowed rather than with elapsed time, the ' ...
    'same mechanism as in tests/checkUS06CoulombCountCapacityMismatch.m. ' ...
    'The EKF is not immune to this either: it uses the same wrong Q_As in ' ...
    'its own prediction step. And just like Experiment C, CC ends up more ' ...
    'accurate here. The 2.9-vs-2.997 Ah mismatch alone is a fairly small ' ...
    'error, while the EKF''s error is still mostly the ~72 mV HPPC-ECM ' ...
    'voltage-fit bias carried over from Experiment A, with this capacity ' ...
    'error only adding a bit on top. CC wins this round, and not because ' ...
    'the EKF failed to use its voltage measurement correctly, that ' ...
    'correction was never going to fix an error baked into the prediction ' ...
    'step itself.\n']);

%% Save results and figure
if ~exist(fullfile(projectRoot, 'results'), 'dir')
    mkdir(fullfile(projectRoot, 'results'));
end
save(fullfile(projectRoot, 'results', 'expE.mat'), ...
    'SOC_ref', 'SOC_CC', 'ekf', 'm_ekf', 'm_cc', 'coverage3Sigma', ...
    'Q_As_measured', 'Q_As_nominal');

figure('Position', [100, 100, 900, 600], 'Color', 'w');
try, set(gcf, 'Theme', 'light'); end  %#ok<TRYNC>  % keep export light even if the MATLAB desktop is in dark mode
subplot(2, 1, 1);
plot(data.t, SOC_ref, 'k-', 'LineWidth', 1.5); hold on;
plot(data.t, ekf.SOC, 'b-');
plot(data.t, SOC_CC, 'r--');
xlabel('Time (s)'); ylabel('SOC (-)');
legend('Reference (measured Ah)', 'EKF (nominal Ah)', 'Coulomb counting (nominal Ah)', 'Location', 'best');
grid on; title('Experiment E - capacity mismatch (documented nominal vs. measured Ah)');

subplot(2, 1, 2);
plot(data.t, err_ekf, 'b-'); hold on;
plot(data.t, err_cc, 'r--');
yline(0, 'k:');
xlabel('Time (s)'); ylabel('SOC error (-)');
legend('EKF error', 'CC error', 'Location', 'best');
grid on; title('CC is more accurate here: EKF error still dominated by Exp A''s model bias');

if ~exist(fullfile(projectRoot, 'figures'), 'dir')
    mkdir(fullfile(projectRoot, 'figures'));
end
exportgraphics(gcf, fullfile(projectRoot, 'figures', 'expE_capacity.png'), 'Resolution', 300, 'BackgroundColor', 'white');
