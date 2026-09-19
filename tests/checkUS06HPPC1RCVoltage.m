%% Fixed HPPC-candidate US06 voltage evidence
% The single Gate-3 comparison for the locked HPPC candidate. Doesn't call
% estimateHPPC1RC, doesn't re-identify or tune any parameter, doesn't touch
% ecmParams, and doesn't use voltage to correct the reference SOC.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

c20File = fullfile(projectRoot, 'data', 'raw', ...
    '05-08-17_13.26 C20 OCV Test_C20_25dC.mat');
us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');

[c20, ~] = loadDataset(string(c20File));
[us06, us06Report] = loadDataset(string(us06File));
[socTarget, voltageTarget] = approvedOCVTarget(c20);
ocvModel = buildOCV(socTarget, voltageTarget);

% Hardcoded straight from the locked one-pass 64-window HPPC run. These are
% validation inputs only, not the parameter set being adopted.
hppcParams = struct;
hppcParams.R0 = 0.0312244650416;
hppcParams.R1 = 0.0113287933401;
hppcParams.tau = 5.23962539484;
Q_Ah = 2.997393193;
Q_As = Q_Ah * 3600;

% Reference trajectory from clean current and the approved initial anchor.
% Voltage doesn't factor into it at all.
SOC_ref = referenceSOC(us06.t, us06.I, 1, Q_As);
assert(all(SOC_ref >= 0 & SOC_ref <= 1), ...
    'The approved US06 reference SOC must remain within the OCV domain.');

[V_model, ~] = ecmSimulate(us06.t, us06.I, SOC_ref, ...
    ocvModel.evaluate, hppcParams);

assert(isnan(V_model(1)), ...
    'The first model-voltage sample must remain undefined.');
comparisonPositions = (2:numel(us06.t)).';
assert(all(isfinite(V_model(comparisonPositions))) && ...
    all(isfinite(us06.V(comparisonPositions))), ...
    'Every aligned comparison sample after the first must be finite.');

voltageError_V = V_model(comparisonPositions) - us06.V(comparisonPositions);
rmse_V = sqrt(mean(voltageError_V .^ 2));
mae_V = mean(abs(voltageError_V));
bias_V = mean(voltageError_V);
maxAbsError_V = max(abs(voltageError_V));

fprintf('\nFixed HPPC-candidate US06 voltage evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  Fixed HPPC candidate: R0 = %.13g ohm, R1 = %.13g ohm, ' ...
    'tau = %.13g s.\n'], hppcParams.R0, hppcParams.R1, hppcParams.tau);
fprintf('  Reference capacity: Q_Ah = %.9f Ah; Q_As = %.9f A*s.\n', ...
    Q_Ah, Q_As);
fprintf(['  V_model(1) = NaN is excluded; aligned samples are positions ' ...
    '2:%d.\n'], numel(us06.t));
fprintf('  Error definition: V_model - V_measured. All metrics are sample-wise.\n');
fprintf('  RMSE: %.9f V\n', rmse_V);
fprintf('  MAE: %.9f V\n', mae_V);
fprintf('  Bias: %+.9f V\n', bias_V);
fprintf('  Maximum absolute error: %.9f V\n', maxAbsError_V);
fprintf(['  This is fixed-parameter evidence only: no adequacy conclusion, tuning, ' ...
    'or voltage-derived SOC correction is made.\n']);

function [socTarget, voltageTarget] = approvedOCVTarget(c20)
% Same approved C/20 target as before, rebuilt here; US06 voltage isn't used.
dischargePositions = find(c20.sourceRow >= 7 & c20.sourceRow <= 1247);
chargePositions = find(c20.sourceRow >= 1309 & c20.sourceRow <= 2391);
dischargeRestPosition = find(c20.sourceRow == 6);
emptyRestPosition = find(c20.sourceRow == 1307);
fullRestPosition = find(c20.sourceRow == 2451);

dischargeDt_s = diff(c20.t([dischargeRestPosition; dischargePositions]));
chargeDt_s = diff(c20.t([emptyRestPosition; chargePositions]));
dischargeCapacity_Ah = ...
    sum(c20.I(dischargePositions) .* dischargeDt_s) / 3600;
chargeCapacity_Ah = ...
    -sum(c20.I(chargePositions) .* chargeDt_s) / 3600;

socDischarge = 1 - [0; cumsum( ...
    c20.I(dischargePositions) .* dischargeDt_s) / 3600] / dischargeCapacity_Ah;
socCharge = [0; cumsum( ...
    -c20.I(chargePositions) .* chargeDt_s) / 3600] / chargeCapacity_Ah;
socDischarge(1) = 1;
socDischarge(end) = 0;
socCharge(1) = 0;
socCharge(end) = 1;

socTarget = linspace(0, 1, 101)';
voltageDischargeGrid = interp1(flipud(socDischarge), ...
    flipud(c20.V([dischargeRestPosition; dischargePositions])), socTarget, ...
    'linear');
voltageChargeGrid = interp1(socCharge, ...
    c20.V([emptyRestPosition; chargePositions]), socTarget, 'linear');

voltageTarget = zeros(size(socTarget));
voltageTarget(1) = c20.V(emptyRestPosition);
voltageTarget(2:100) = ...
    (voltageDischargeGrid(2:100) + voltageChargeGrid(2:100)) / 2;
voltageTarget(end) = mean([c20.V(dischargeRestPosition), ...
    c20.V(fullRestPosition)]);
end
