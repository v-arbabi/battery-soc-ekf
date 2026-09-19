function tests = testLoadDataset
%TESTLOADDATASET Verifies both approved dataset files load correctly.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.projectRoot = projectRoot;
testCase.TestData.rawDataDir = fullfile(projectRoot, 'data', 'raw');
addpath(fullfile(projectRoot, 'src'));
end

function testC20SelectionSignDuplicatesAndTime(testCase)
fileName = "05-08-17_13.26 C20 OCV Test_C20_25dC.mat";
filePath = fullfile(testCase.TestData.rawDataDir, fileName);
raw = load(filePath, 'meas');

[data, report] = loadDataset(filePath);

expectedSelectedRows = (1:2452)';
expectedRemovedRows = [1308; 2452];
expectedKeptRows = setdiff(expectedSelectedRows, expectedRemovedRows, 'stable');

verifyEqual(testCase, report.selectedSourceRows, expectedSelectedRows);
verifyEqual(testCase, report.excludedSourceRows, 2453);
verifyEqual(testCase, report.removedDuplicateSourceRows, expectedRemovedRows);
verifyEqual(testCase, data.sourceRow, expectedKeptRows);
verifyEqual(testCase, data.t, raw.meas.Time(expectedKeptRows));
verifyEqual(testCase, data.I, -raw.meas.Current(expectedKeptRows));
verifyEqual(testCase, data.V, raw.meas.Voltage(expectedKeptRows));
verifyEqual(testCase, data.T, raw.meas.Battery_Temp_degC(expectedKeptRows));
verifyEqual(testCase, diff(data.t), diff(raw.meas.Time(expectedKeptRows)));
verifyTrue(testCase, all(diff(data.t) > 0));
verifyFalse(testCase, report.interpolationApplied);
verifyEqual(testCase, report.currentNormalisation, ...
    "I = -meas.Current (applied once)");
end

function testUS06SignDuplicatesAndNonuniformTime(testCase)
fileName = "03-20-17_01.43 25degC_US06_Pan18650PF.mat";
filePath = fullfile(testCase.TestData.rawDataDir, fileName);
raw = load(filePath, 'meas');

[data, report] = loadDataset(filePath);

expectedSelectedRows = (1:48061)';
expectedRemovedRows = 48061;
expectedKeptRows = (1:48060)';
expectedDt = diff(raw.meas.Time(expectedKeptRows));

verifyEqual(testCase, report.selectedSourceRows, expectedSelectedRows);
verifyEmpty(testCase, report.excludedSourceRows);
verifyEqual(testCase, report.removedDuplicateSourceRows, expectedRemovedRows);
verifyEqual(testCase, data.sourceRow, expectedKeptRows);
verifyEqual(testCase, data.t, raw.meas.Time(expectedKeptRows));
verifyEqual(testCase, data.I, -raw.meas.Current(expectedKeptRows));
verifyEqual(testCase, data.V, raw.meas.Voltage(expectedKeptRows));
verifyEqual(testCase, data.T, raw.meas.Battery_Temp_degC(expectedKeptRows));
verifyEqual(testCase, diff(data.t), expectedDt);
verifyTrue(testCase, any(expectedDt > 1.8));
verifyGreaterThan(testCase, numel(unique(expectedDt)), 1);
verifyTrue(testCase, all(diff(data.t) > 0));
verifyFalse(testCase, report.interpolationApplied);
verifyEqual(testCase, report.currentNormalisation, ...
    "I = -meas.Current (applied once)");
end

function testRejectsDecreasingTimeWithExactIdentifier(testCase)
sourceFile = fullfile(testCase.TestData.rawDataDir, ...
    "05-08-17_13.26 C20 OCV Test_C20_25dC.mat");
loadedData = load(sourceFile, 'meas');
meas = loadedData.meas;
meas.Time(100) = meas.Time(99) - 1;

fixture = matlab.unittest.fixtures.TemporaryFolderFixture;
testCase.applyFixture(fixture);
temporaryFile = fullfile(fixture.Folder, ...
    "05-08-17_13.26 C20 OCV Test_C20_25dC.mat");
save(temporaryFile, 'meas');

verifyError(testCase, @() loadDataset(temporaryFile), ...
    "loadDataset:DecreasingTime");
end

function testRejectsNonIdenticalZeroDtWithExactIdentifier(testCase)
sourceFile = fullfile(testCase.TestData.rawDataDir, ...
    "05-08-17_13.26 C20 OCV Test_C20_25dC.mat");
loadedData = load(sourceFile, 'meas');
meas = loadedData.meas;
meas.Time(100) = meas.Time(99);
meas.Voltage(100) = meas.Voltage(99) + 1e-6;

fixture = matlab.unittest.fixtures.TemporaryFolderFixture;
testCase.applyFixture(fixture);
temporaryFile = fullfile(fixture.Folder, ...
    "05-08-17_13.26 C20 OCV Test_C20_25dC.mat");
save(temporaryFile, 'meas');

verifyError(testCase, @() loadDataset(temporaryFile), ...
    "loadDataset:NonIdenticalZeroDt");
end
