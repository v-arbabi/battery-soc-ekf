%% US06 Coulomb Counting documented-capacity-mismatch evidence
% Current-only check showing the gap between the documented nominal
% capacity and the measured one. This isn't SOH estimation, isn't a claim
% about a degraded cell, and isn't external validation. It just shows what
% that capacity assumption does to SOC_CC.

clearvars;

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));

us06File = fullfile(projectRoot, 'data', 'raw', ...
    '03-20-17_01.43 25degC_US06_Pan18650PF.mat');
[us06, us06Report] = loadDataset(string(us06File));

SOC0 = 1;
measuredQ_Ah = 2.997393193;
nominalQ_Ah = 2.9;
measuredQ_As = measuredQ_Ah * 3600;
nominalQ_As = nominalQ_Ah * 3600;

% Reference SOC comes from clean current and the measured capacity.
SOC_ref = referenceSOC(us06.t, us06.I, SOC0, measuredQ_As);

% Baseline only sees clean time/current, its initial SOC, and the nominal
% capacity from the datasheet.
SOC_CC = coulombCount(us06.t, us06.I, SOC0, nominalQ_As);

comparison = metrics(SOC_CC, SOC_ref);
initialError = comparison.error(1);
finalError = comparison.error(end);

fprintf('\nUS06 Coulomb Counting documented-capacity-mismatch evidence\n');
fprintf('  Retained US06 samples: %d. Removed duplicate source rows: %s.\n', ...
    numel(us06.t), mat2str(us06Report.removedDuplicateSourceRows.'));
fprintf(['  Shared clean inputs: actual numeric time, sign-normalised current, ' ...
    'and SOC0 = %.1f.\n'], SOC0);
fprintf(['  Reference capacity: %.9f Ah (%.9f A*s); baseline assumed capacity: ' ...
    '%.9f Ah (%.9f A*s).\n'], measuredQ_Ah, measuredQ_As, ...
    nominalQ_Ah, nominalQ_As);
fprintf('  Error definition: SOC_CC - SOC_ref. All metrics are sample-wise.\n');
fprintf('  Aligned valid samples: %d\n', comparison.validSampleCount);
fprintf('  RMSE: %.12g\n', comparison.RMSE);
fprintf('  MAE: %.12g\n', comparison.MAE);
fprintf('  Bias: %+.12g\n', comparison.bias);
fprintf('  Maximum absolute error: %.12g\n', comparison.maxAbsError);
fprintf('  Initial error: %+.12g\n', initialError);
fprintf('  Final error: %+.12g\n', finalError);
fprintf(['  This is a documented nominal-versus-measured capacity-assumption ' ...
    'mismatch only; it is not SOH estimation or a new degraded-plant claim.\n']);
