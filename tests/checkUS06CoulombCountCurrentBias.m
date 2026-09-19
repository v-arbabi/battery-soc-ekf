%% US06 Coulomb Counting fixed-current-bias evidence
% Adds a fixed illustrative current bias and watches SOC_CC drift away
% from SOC_ref. Not a real sensor characterisation, not external
% validation, and not a comparison against the EKF.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[us06, us06Report] = loadDataset(string(us06File));

SOC0 = 1;
Q_Ah = 2.997393193;
Q_As = Q_Ah * 3600;
currentBias_A = 0.02;
I_clean = us06.I;
I_biased = I_clean + currentBias_A;

% Reference SOC still comes from the clean current, untouched.
SOC_ref = referenceSOC(us06.t, I_clean, SOC0, Q_As);

% Baseline only takes actual time, the biased current, its initial SOC,
% and capacity. No reference, voltage, OCV, or ECM value goes into it.
SOC_CC = coulombCount(us06.t, I_biased, SOC0, Q_As);

comparison = metrics(SOC_CC, SOC_ref);
initialError = comparison.error(1);
finalError = comparison.error(end);

fprintf('\nUS06 Coulomb Counting fixed-current-bias evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  Illustrative experiment assumption: I_biased = I_clean + %.3f A ' ...
    'at every retained sample.\n'], currentBias_A);
fprintf('  Shared SOC0 = %.1f and Q_As = %.9f A*s.\n', SOC0, Q_As);
fprintf('  Error definition: SOC_CC - SOC_ref. All metrics are sample-wise.\n');
fprintf('  Aligned valid samples: %d\n', comparison.validSampleCount);
fprintf('  RMSE: %.12g\n', comparison.RMSE);
fprintf('  MAE: %.12g\n', comparison.MAE);
fprintf('  Bias: %+.12g\n', comparison.bias);
fprintf('  Maximum absolute error: %.12g\n', comparison.maxAbsError);
fprintf('  Initial error: %+.12g\n', initialError);
fprintf('  Final error: %+.12g\n', finalError);
fprintf(['  A positive project-sign current bias produces cumulative negative ' ...
    'SOC_CC - SOC_ref drift under positive-discharge integration.\n']);
fprintf(['  The +0.02 A value is an illustrative experiment assumption, not a ' ...
    'measured sensor specification or external-accuracy claim.\n']);
