%% US06 matched-integration Coulomb Counting consistency evidence
% Read-only check: builds SOC_CC and SOC_ref from the same clean
% current/time inputs and compares them. Just an internal consistency
% check, not external validation and not a claim about estimator accuracy.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[us06, us06Report] = loadDataset(string(us06File));

SOC0 = 1;
Q_Ah = 2.997393193;
Q_As = Q_Ah * 3600;

% coulombCount only gets the retained time/current, the initial SOC, and
% the capacity here. Reference SOC never goes near this call.
SOC_CC = coulombCount(us06.t, us06.I, SOC0, Q_As);

% Reference trajectory built separately, same clean inputs and the same
% right-endpoint integration rule.
SOC_ref = referenceSOC(us06.t, us06.I, SOC0, Q_As);

comparison = metrics(SOC_CC, SOC_ref);
assert(comparison.validSampleCount == numel(us06.t), ...
    'Every retained US06 position must enter the matched-integration comparison.');

fprintf('\nUS06 Coulomb Counting matched-integration consistency evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  Both integrations use clean, sign-normalised US06 current, actual ' ...
    'numeric time, SOC0 = %.1f, and Q_As = %.9f A*s.\n'], SOC0, Q_As);
fprintf('  SOC_CC final value: %.12f\n', SOC_CC(end));
fprintf('  SOC_ref final value: %.12f\n', SOC_ref(end));
fprintf('  Error definition: SOC_CC - SOC_ref. All metrics are sample-wise.\n');
fprintf('  Aligned valid samples: %d\n', comparison.validSampleCount);
fprintf('  RMSE: %.12g\n', comparison.RMSE);
fprintf('  MAE: %.12g\n', comparison.MAE);
fprintf('  Bias: %+.12g\n', comparison.bias);
fprintf('  Maximum absolute error: %.12g\n', comparison.maxAbsError);
fprintf(['  This is matched-integration consistency evidence only; it does not ' ...
    'establish external accuracy or estimator performance.\n']);
