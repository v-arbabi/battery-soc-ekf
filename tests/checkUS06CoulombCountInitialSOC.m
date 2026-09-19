%% US06 Coulomb Counting incorrect-initial-SOC evidence
% Starts the baseline from a wrong initial SOC and tracks how the error
% behaves. Says nothing about external accuracy, recovery, or how an
% estimator would actually perform.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[us06, us06Report] = loadDataset(string(us06File));

referenceSOC0 = 1;
baselineSOC0 = 0.8;
Q_Ah = 2.997393193;
Q_As = Q_Ah * 3600;

% Baseline goes first: retained clean time/current, its own (wrong) initial
% SOC, and the approved capacity.
SOC_CC = coulombCount(us06.t, us06.I, baselineSOC0, Q_As);

% Reference SOC built separately, and coulombCount never sees it. The two
% trajectories only come together down in metrics() below.
SOC_ref = referenceSOC(us06.t, us06.I, referenceSOC0, Q_As);

comparison = metrics(SOC_CC, SOC_ref);
initialError = comparison.error(1);
finalError = comparison.error(end);
maximumOffsetDeviation = max(abs(comparison.error - initialError));
fixedOffsetPreserved = maximumOffsetDeviation <= 1e-12;

fprintf('\nUS06 Coulomb Counting incorrect-initial-SOC evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  Shared clean inputs: actual numeric time, sign-normalised current, ' ...
    'and Q_As = %.9f A*s.\n'], Q_As);
fprintf('  Reference initial SOC: %.1f; baseline initial SOC: %.1f.\n', ...
    referenceSOC0, baselineSOC0);
fprintf('  Error definition: SOC_CC - SOC_ref. All metrics are sample-wise.\n');
fprintf('  Aligned valid samples: %d\n', comparison.validSampleCount);
fprintf('  RMSE: %.12g\n', comparison.RMSE);
fprintf('  MAE: %.12g\n', comparison.MAE);
fprintf('  Bias: %+.12g\n', comparison.bias);
fprintf('  Maximum absolute error: %.12g\n', comparison.maxAbsError);
fprintf('  Initial error: %+.12g\n', initialError);
fprintf('  Final error: %+.12g\n', finalError);
fprintf(['  Fixed initial offset preserved within 1e-12: %s ' ...
    '(maximum deviation %.12g).\n'], string(fixedOffsetPreserved), ...
    maximumOffsetDeviation);
fprintf(['  This isolates the declared initial-SOC mismatch under otherwise ' ...
    'matched integration; it does not establish external accuracy, recovery, ' ...
    'or estimator performance.\n']);
