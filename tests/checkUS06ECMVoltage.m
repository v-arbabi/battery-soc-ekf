%% Focused US06 forward-model voltage evidence
% Runs the approved loader, current-integrated reference SOC, degree-9
% C/20 OCV fit, and the fixed provisional ECM parameters through the
% right-endpoint simulator. Measured US06 voltage only enters afterward,
% to compute error stats. No parameter tuning or SOC correction happens here.

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
params = ecmParams();

% Protocol-supported anchor, clean-current integration. Voltage plays no
% part in this calculation.
SOC_ref = referenceSOC(us06.t, us06.I, 1, params.Q_As);
assert(all(SOC_ref >= 0 & SOC_ref <= 1), ...
    'The approved US06 reference SOC must remain within the OCV domain.');

[V_model, V1] = ecmSimulate(us06.t, us06.I, SOC_ref, ...
    ocvModel.evaluate, params);

assert(isnan(V_model(1)), ...
    'The first model-voltage sample must remain undefined.');
validPositions = find(isfinite(V_model) & isfinite(us06.V));
assert(isequal(validPositions, (2:numel(us06.t))'), ...
    'Only aligned samples 2:end may enter US06 voltage statistics.');

voltageError_V = V_model(validPositions) - us06.V(validPositions);
rmse_V = sqrt(mean(voltageError_V .^ 2));
mae_V = mean(abs(voltageError_V));
bias_V = mean(voltageError_V);
maxAbsError_V = max(abs(voltageError_V));

alignedCurrent_A = us06.I(validPositions);
alignedSOC = SOC_ref(validPositions);

fprintf('\nUS06 forward-model voltage evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  First model-voltage sample: V_model(1) = NaN; it is excluded ' ...
    'because no preceding retained interval exists.\n']);
fprintf('  Aligned finite comparison samples: %d (source-model positions 2:%d).\n', ...
    numel(validPositions), numel(us06.t));
fprintf(['  Inputs held fixed: R0 = %.9f ohm, R1 = %.9f ohm, tau = %.9f s, ' ...
    'and V1(1) = %.1f V.\n'], ...
    params.R0, params.R1, params.tau, V1(1));
fprintf(['  SOC_ref uses clean current, Q_As = %.7f A*s, and the explicit ' ...
    'protocol-supported SOC_ref(1) = %.1f; it is not voltage-corrected.\n'], ...
    params.Q_As, SOC_ref(1));
fprintf('\nVoltage error definition: V_model - V_measured\n');
fprintf('  RMSE: %.9f V\n', rmse_V);
fprintf('  MAE: %.9f V\n', mae_V);
fprintf('  Bias: %+.9f V\n', bias_V);
fprintf('  Maximum absolute error: %.9f V\n', maxAbsError_V);
fprintf(['  All aggregate metrics in this check are sample-wise; no metric is ' ...
    'duration-weighted.\n']);

fprintf('\nDescriptive residual groups versus signed current (sample-wise)\n');
currentEdges_A = [-Inf, -0.5, -0.05, 0.05, 0.5, Inf];
currentLabels = [ ...
    "I < -0.50 A"
    "-0.50 <= I < -0.05 A"
    "-0.05 <= I < 0.05 A"
    "0.05 <= I < 0.50 A"
    "I >= 0.50 A"
];
printResidualGroups(alignedCurrent_A, voltageError_V, currentEdges_A, ...
    currentLabels);

fprintf('\nDescriptive residual groups versus SOC_ref (sample-wise)\n');
socEdges = [0, 0.2, 0.4, 0.6, 0.8, 1];
socLabels = [ ...
    "0.00 <= SOC_ref < 0.20"
    "0.20 <= SOC_ref < 0.40"
    "0.40 <= SOC_ref < 0.60"
    "0.60 <= SOC_ref < 0.80"
    "0.80 <= SOC_ref <= 1.00"
];
printResidualGroups(alignedSOC, voltageError_V, socEdges, socLabels);
fprintf(['\nThese are descriptive evidence for the fixed provisional model only; ' ...
    'no physical-cause attribution, parameter adequacy conclusion, tuning, ' ...
    'or SOC correction is made.\n']);

function [socTarget, voltageTarget] = approvedOCVTarget(c20)
% Rebuilds the already-approved C/20 target; doesn't touch US06 voltage.
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

function printResidualGroups(groupValue, voltageError_V, edges, labels)
for groupIndex = 1:numel(labels)
    if groupIndex == numel(labels)
        inGroup = groupValue >= edges(groupIndex) & ...
            groupValue <= edges(groupIndex + 1);
    else
        inGroup = groupValue >= edges(groupIndex) & ...
            groupValue < edges(groupIndex + 1);
    end

    groupError_V = voltageError_V(inGroup);
    if isempty(groupError_V)
        fprintf('  %-28s n = %6d; no samples.\n', labels(groupIndex), 0);
    else
        fprintf(['  %-28s n = %6d; mean = %+.9f V; RMSE = %.9f V; ' ...
            'MAE = %.9f V; max |error| = %.9f V.\n'], ...
            labels(groupIndex), numel(groupError_V), mean(groupError_V), ...
            sqrt(mean(groupError_V .^ 2)), mean(abs(groupError_V)), ...
            max(abs(groupError_V)));
    end
end
end
