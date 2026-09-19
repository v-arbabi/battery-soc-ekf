function [socTarget, voltageTarget] = buildOCVTarget(projectRoot)
%BUILDOCVTARGET Construct the approved 101-point OCV-SOC target.
%   [SOCTARGET, VOLTAGETARGET] = BUILDOCVTARGET(PROJECTROOT) rebuilds the
%   branch-mean OCV target from the C/20 characterisation file (see
%   docs/DECISIONS.md, "OCV(SOC) model"). This logic previously lived as a
%   local function inside tests/testBuildOCV.m and was pulled out here so
%   other code can call it too. Nothing about the anchors, boundaries, or
%   endpoint policy changed in the move.
%
%   Anchors, unchanged from the approved decision:
%     SOC 0: rested source row 1307 only.
%     SOC 1: mean of rested source rows 6 and 2451.
%     Interior (SOC 0.01-0.99): arithmetic mean of the separately
%       normalised discharge and charge branches. This is low-rate
%       charge/discharge voltage separation, not equilibrium hysteresis.
%     Loaded source rows 1247 and 2391 stay out of the endpoint anchors.

arguments
    projectRoot (1, 1) string
end

c20File = fullfile(projectRoot, 'data', 'raw', ...
    '05-08-17_13.26 C20 OCV Test_C20_25dC.mat');
[data, ~] = loadDataset(string(c20File));

dischargePositions = find(data.sourceRow >= 7 & data.sourceRow <= 1247);
chargePositions = find(data.sourceRow >= 1309 & data.sourceRow <= 2391);
dischargeRestPosition = find(data.sourceRow == 6);
emptyRestPosition = find(data.sourceRow == 1307);
fullRestPosition = find(data.sourceRow == 2451);

dischargeDt = diff(data.t([dischargeRestPosition; dischargePositions]));
chargeDt = diff(data.t([emptyRestPosition; chargePositions]));
dischargeCapacityAh = ...
    sum(data.I(dischargePositions) .* dischargeDt) / 3600;
chargeCapacityAh = ...
    -sum(data.I(chargePositions) .* chargeDt) / 3600;

socDischarge = 1 - [0; cumsum( ...
    data.I(dischargePositions) .* dischargeDt) / 3600] / dischargeCapacityAh;
socCharge = [0; cumsum( ...
    -data.I(chargePositions) .* chargeDt) / 3600] / chargeCapacityAh;
socDischarge(1) = 1;
socDischarge(end) = 0;
socCharge(1) = 0;
socCharge(end) = 1;

voltageDischarge = data.V([dischargeRestPosition; dischargePositions]);
voltageCharge = data.V([emptyRestPosition; chargePositions]);
socTarget = linspace(0, 1, 101)';
voltageDischargeGrid = interp1(flipud(socDischarge), ...
    flipud(voltageDischarge), socTarget, 'linear');
voltageChargeGrid = interp1(socCharge, voltageCharge, socTarget, 'linear');

voltageTarget = zeros(size(socTarget));
voltageTarget(1) = data.V(emptyRestPosition);
voltageTarget(2:100) = ...
    (voltageDischargeGrid(2:100) + voltageChargeGrid(2:100)) / 2;
voltageTarget(101) = mean([data.V(dischargeRestPosition), ...
    data.V(fullRestPosition)]);
end
