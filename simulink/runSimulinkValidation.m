%RUNSIMULINKVALIDATION Build (if needed) and validate ecm_ekf.slx against
%the MATLAB reference src/ekfEstimate.m, under the same Experiment A
%conditions (nominal, HPPC parameters, clean US06 current and voltage).
%
%   This script always calls buildModel() first, so it tests whatever
%   buildModel.m currently says, never a stale .slx left over from an
%   earlier run. Confirmed on 2026-09-18:
%   max |SOC_simulink - SOC_ekfEstimate| = 0.000e+00.

clear; clc;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
addpath(fullfile(projectRoot, 'simulink'));

% Rebuild from buildModel.m every time rather than trusting whatever .slx
% is sitting on disk. buildModel.m has already been corrected more than
% once during testing, and a stale .slx could quietly reproduce a bug
% that's already been fixed.
buildModel();

%% Load the same data and OCV model as experiments/expA_nominal.m
us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[data, ~] = loadDataset(string(us06File));

[socTarget, voltageTarget] = buildOCVTarget(string(projectRoot));
ocvModel = buildOCV(socTarget, voltageTarget);

paramsHPPC = ecmParamsHPPC();

%% Reference run: the approved MATLAB EKF, Experiment A tuning
tuning = struct;
tuning.P0_SOC = (0.01)^2;
tuning.P0_V1  = (0.01)^2;
tuning.q_SOC  = 1e-10;
tuning.q_V1   = 1e-6;
tuning.R_V    = (3e-3)^2;

ekfRef = ekfEstimate(data.t, data.I, data.V, ...
    ocvModel.evaluate, ocvModel.evaluateDerivative, ...
    paramsHPPC, 1.0, paramsHPPC.Q_As, tuning);

%% Prepare base-workspace signals for the Simulink model's From Workspace
%  blocks. Time base is the plain sample index (0, 1, ..., N-1), not real
%  seconds - see the correction note in buildModel.m for the reason: it's
%  what forces Simulink to take exactly one fixed-step tick per data row
%  instead of trying (and failing) to hit irregular real-time seconds.
%  Real elapsed time still gets carried as data, in dt_ws's second column;
%  dt_ws(1,2) = 0 is the first-sample flag the "EKF Step" MATLAB Function
%  block checks for (see buildModel.m).
N = numel(data.t);
tick = (0:N-1)';
dt_col = [0; diff(data.t)];
dt_ws = [tick, dt_col];
I_ws  = [tick, data.I];
V_ws  = [tick, data.V];

ocvCoeffs_ws  = ocvModel.coefficients;           %#ok<NASGU> % degree-9, 10 coeffs
docvCoeffs_ws = ocvModel.derivativeCoefficients; %#ok<NASGU>

assignin('base', 'dt_ws', dt_ws);
assignin('base', 'I_ws', I_ws);
assignin('base', 'V_ws', V_ws);
assignin('base', 'ocvCoeffs_ws', ocvCoeffs_ws);
assignin('base', 'docvCoeffs_ws', docvCoeffs_ws);

%% Run the Simulink model for exactly N-1 ticks (N samples, tick 0 to N-1)
simOut = sim('ecm_ekf', 'StopTime', num2str(N - 1));

% Newer MATLAB versions return To Workspace outputs as fields on the
% Simulink.SimulationOutput object (out.SOC_simulink, matching the
% SOC_out block name in the diagram) rather than as base-workspace
% variables. Try that first, then fall back to the base workspace for
% older versions that still populate it directly.
if isprop(simOut, 'SOC_simulink') || (ismethod(simOut, 'get') && any(strcmp(simOut.who, 'SOC_simulink')))
    SOC_simulink = simOut.SOC_simulink;
elseif evalin('base', 'exist(''SOC_simulink'', ''var'')')
    SOC_simulink = evalin('base', 'SOC_simulink');
else
    error('runSimulinkValidation:CannotFindOutput', ...
        ['Could not find SOC_simulink in simOut or the base workspace. ' ...
        'Run "simOut" alone at the command line right after this line and ' ...
        'look at what fields/properties it lists, then report that back.']);
end

%% Compare
assert(numel(SOC_simulink) == N, 'runSimulinkValidation:SizeMismatch', ...
    'Simulink logged %d samples, expected %d. Check From Workspace timing.', ...
    numel(SOC_simulink), N);

diffSOC = SOC_simulink(:) - ekfRef.SOC;
maxAbsDiff = max(abs(diffSOC));

fprintf('\nSimulink vs. MATLAB (src/ekfEstimate.m) EKF SOC comparison, Experiment A conditions\n');
fprintf('  max |SOC_simulink - SOC_ekfEstimate| = %.3e\n', maxAbsDiff);

tol = 1e-9;
if maxAbsDiff <= tol
    fprintf('  PASS: Simulink model matches the MATLAB reference to %.0e.\n', tol);
else
    fprintf(['  FAIL: difference exceeds %.0e. Do not treat the Simulink model as ' ...
        'validated - report the max difference and where it first appears ' ...
        '(min(find(abs(diffSOC) > tol))) so the block logic can be corrected.\n'], tol);
end
