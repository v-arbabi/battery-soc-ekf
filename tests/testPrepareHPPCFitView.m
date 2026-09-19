function tests = testPrepareHPPCFitView
%TESTPREPAREHPPCFITVIEW Strict-time HPPC provenance checks.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'src'));
testCase.TestData.projectRoot = projectRoot;
end

function testExactAndNonidenticalZeroDtPolicy(testCase)
meas = syntheticMeas();
originalMeas = meas;

fitView = prepareHPPCFitView(meas);

verifyEqual(testCase, meas, originalMeas);
verifyEqual(testCase, fitView.sourceRow, [1; 2; 6]);
verifyEqual(testCase, fitView.report.retainedSourceRows, [1; 2; 6]);
verifyEqual(testCase, fitView.report.exactDuplicateRemovedSourceRows, 3);
verifyEqual(testCase, fitView.report.nonidenticalExcludedSourceRows, [4; 5]);
verifyEqual(testCase, fitView.t, [0; 1; 4]);
verifyEqual(testCase, diff(fitView.t), [1; 3]);
verifyTrue(testCase, all(diff(fitView.t) > 0));
verifyEqual(testCase, fitView.rawCurrent, [0; -2; 0]);
verifyEqual(testCase, fitView.I, [0; 2; 0]);
verifyEqual(testCase, fitView.V, [4.2; 4.0; 3.9]);
verifyEqual(testCase, fitView.report.currentNormalisation, ...
    "I = -meas.Current (applied once)");
verifyFalse(testCase, fitView.report.interpolationApplied);
verifyFalse(testCase, fitView.report.timeRepairApplied);
end

function testNoZeroDtRowsAreRetainedWithoutTimeRepair(testCase)
meas = syntheticMeas();
meas = measSubset(meas, [1; 2; 6]);

fitView = prepareHPPCFitView(meas);

verifyEqual(testCase, fitView.sourceRow, [1; 2; 3]);
verifyEqual(testCase, fitView.t, [0; 1; 4]);
verifyEqual(testCase, diff(fitView.t), [1; 3]);
verifyEmpty(testCase, fitView.report.exactDuplicateRemovedSourceRows);
verifyEmpty(testCase, fitView.report.nonidenticalExcludedSourceRows);
end

function testFullHPPCProvenanceAndStrictTime(testCase)
hppcFile = fullfile(testCase.TestData.projectRoot, 'data', 'raw', ...
    '03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat');
loaded = load(hppcFile, 'meas');

fitView = prepareHPPCFitView(loaded.meas);

verifyEqual(testCase, fitView.report.rawSampleCount, 102800);
verifyEqual(testCase, numel(fitView.report.exactDuplicateRemovedSourceRows), 133);
verifyEqual(testCase, numel(fitView.report.nonidenticalExcludedSourceRows), 40);
verifyEqual(testCase, numel(fitView.sourceRow), 102627);
verifyTrue(testCase, all(diff(fitView.t) > 0));
verifyTrue(testCase, all(fitView.I == -fitView.rawCurrent));
verifyFalse(testCase, any(fitView.sourceRow == 202));
verifyTrue(testCase, any(fitView.sourceRow == 201));
verifyFalse(testCase, any(ismember(fitView.sourceRow, [4488; 4489])));
verifyEmpty(testCase, intersect(fitView.report.exactDuplicateRemovedSourceRows, ...
    fitView.report.nonidenticalExcludedSourceRows));
end

function testInputValidationErrorIdentifiers(testCase)
meas = syntheticMeas();

verifyError(testCase, @() prepareHPPCFitView(struct()), ...
    "prepareHPPCFitView:InvalidMeas");
missingField = rmfield(meas, 'Power');
verifyError(testCase, @() prepareHPPCFitView(missingField), ...
    "prepareHPPCFitView:MissingFields");
inconsistentLength = meas;
inconsistentLength.Voltage = inconsistentLength.Voltage(1:end-1);
verifyError(testCase, @() prepareHPPCFitView(inconsistentLength), ...
    "prepareHPPCFitView:InconsistentFieldLength");
nonfiniteSignal = meas;
nonfiniteSignal.Current(1) = NaN;
verifyError(testCase, @() prepareHPPCFitView(nonfiniteSignal), ...
    "prepareHPPCFitView:NonFiniteSignal");
decreasingTime = measSubset(meas, [1; 2; 6]);
decreasingTime.Time = [0; 2; 1];
verifyError(testCase, @() prepareHPPCFitView(decreasingTime), ...
    "prepareHPPCFitView:DecreasingTime");
end

function testAmbiguousOverlappingZeroDtRunIsRejected(testCase)
meas = syntheticMeas();
% Rows 2-3 are exact duplicates; rows 3-4 share the same timestamp but
% aren't identical. That puts row 3 in both zero-dt buckets at once, which
% is exactly the ambiguous case that should get rejected outright.
meas.Time = [0; 1; 1; 1; 2.7; 4];

verifyError(testCase, @() prepareHPPCFitView(meas), ...
    "prepareHPPCFitView:AmbiguousZeroDtRun");
end

function meas = syntheticMeas()
sampleCount = 6;
meas = struct;
meas.TimeStamp = { ...
    '01/01/2026 12:00:00 AM'; '01/01/2026 12:00:01 AM'; ...
    '01/01/2026 12:00:01 AM'; '01/01/2026 12:00:02 AM'; ...
    '01/01/2026 12:00:02 AM'; '01/01/2026 12:00:04 AM'};
meas.Voltage = [4.2; 4.0; 4.0; 3.5; 3.8; 3.9];
meas.Current = [0; -2; -2; -4; 0; 0];
meas.Ah = [0; -0.01; -0.01; -0.02; -0.02; -0.02];
meas.Wh = [0; -0.04; -0.04; -0.07; -0.07; -0.07];
meas.Power = meas.Voltage .* meas.Current;
meas.Battery_Temp_degC = 25 * ones(sampleCount, 1);
meas.Time = [0; 1; 1; 2.7; 2.7; 4];
meas.Chamber_Temp_degC = 25 * ones(sampleCount, 1);
end

function subset = measSubset(meas, rows)
fields = fieldnames(meas);
subset = struct;
for fieldIndex = 1:numel(fields)
    fieldName = fields{fieldIndex};
    subset.(fieldName) = meas.(fieldName)(rows);
end
end
