function [data, report] = loadDataset(filePath)
%LOADDATASET Load one of the two approved Panasonic 18650PF dataset files.
%   [DATA, REPORT] = LOADDATASET(FILEPATH) loads the scalar meas struct,
%   applies the approved source-row boundary, drops the second row of any
%   exact zero-dt duplicate pair, and normalises current once (positive =
%   discharge). Real time intervals stay as-is, no interpolation, no
%   swapping in a constant dt.
%
%   DATA fields:
%     t         time, seconds
%     I         current, amps, positive for discharge
%     V         terminal voltage, volts
%     T         battery case temperature, degrees C
%     sourceRow original one-based row number in the MAT file
%
%   REPORT records which source rows were kept and which duplicates got
%   removed.

arguments
    filePath (1, 1) string
end

c20FileName = "05-08-17_13.26 C20 OCV Test_C20_25dC.mat";
us06FileName = "03-20-17_01.43 25degC_US06_Pan18650PF.mat";
requiredFields = [
    "TimeStamp"
    "Voltage"
    "Current"
    "Ah"
    "Wh"
    "Power"
    "Battery_Temp_degC"
    "Time"
    "Chamber_Temp_degC"
];

assert(isfile(filePath), "loadDataset:FileNotFound", ...
    "Dataset file not found: %s", filePath);

loadedData = load(filePath, 'meas');
assert(isfield(loadedData, 'meas') && isscalar(loadedData.meas), ...
    "loadDataset:InvalidMeas", ...
    "Expected one scalar struct named meas in %s.", filePath);
meas = loadedData.meas;

missingFields = setdiff(requiredFields, string(fieldnames(meas)));
assert(isempty(missingFields), "loadDataset:MissingFields", ...
    "Missing required fields: %s", join(missingFields, ", "));

sampleCount = numel(meas.Time);
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields(fieldIndex);
    assert(isvector(meas.(fieldName)) && ...
        numel(meas.(fieldName)) == sampleCount, ...
        "loadDataset:InconsistentFieldLength", ...
        "Field %s must be a vector with %d elements.", fieldName, sampleCount);
end

[~, fileStem, fileExtension] = fileparts(filePath);
fileName = string(fileStem) + string(fileExtension);
if fileName == c20FileName
    assert(sampleCount == 2453, "loadDataset:UnexpectedC20Length", ...
        "The approved C/20 source file must contain 2453 rows.");
    selectedSourceRows = (1:2452)';
elseif fileName == us06FileName
    selectedSourceRows = (1:sampleCount)';
else
    error("loadDataset:UnsupportedFile", ...
        "File is not one of the two approved datasets: %s", fileName);
end

selectedTime = meas.Time(selectedSourceRows);
assert(isnumeric(selectedTime) && all(isfinite(selectedTime)), ...
    "loadDataset:InvalidTime", "Time must contain finite numeric values.");

dt = diff(selectedTime);
assert(~any(dt < 0), "loadDataset:DecreasingTime", ...
    "Numeric Time decreases within the approved source rows.");

duplicatePositions = find(dt == 0) + 1;
removeMask = false(size(selectedSourceRows));
removedDuplicateSourceRows = zeros(numel(duplicatePositions), 1);
for duplicateIndex = 1:numel(duplicatePositions)
    secondPosition = duplicatePositions(duplicateIndex);
    firstSourceRow = selectedSourceRows(secondPosition - 1);
    secondSourceRow = selectedSourceRows(secondPosition);

    nonIdenticalZeroDtMessage = sprintf( ...
        ['Rows %d and %d have zero dt but differ in at least one ' ...
        'required field.'], firstSourceRow, secondSourceRow);
    assert(rowsAreExactDuplicates(meas, requiredFields, ...
        firstSourceRow, secondSourceRow), ...
        "loadDataset:NonIdenticalZeroDt", nonIdenticalZeroDtMessage);

    removeMask(secondPosition) = true;
    removedDuplicateSourceRows(duplicateIndex) = secondSourceRow;
end

keptSourceRows = selectedSourceRows(~removeMask);
keptTime = meas.Time(keptSourceRows);
assert(all(diff(keptTime) > 0), "loadDataset:InvalidCleanedTime", ...
    "Time must be strictly increasing after exact duplicate removal.");

data = struct;
data.t = keptTime(:);
data.I = -meas.Current(keptSourceRows);
data.V = meas.Voltage(keptSourceRows);
data.T = meas.Battery_Temp_degC(keptSourceRows);
data.sourceRow = keptSourceRows;

data.I = data.I(:);
data.V = data.V(:);
data.T = data.T(:);

assert(all(isfinite(data.I)) && all(isfinite(data.V)) && all(isfinite(data.T)), ...
    "loadDataset:NonFiniteSignal", ...
    "Current, voltage, and battery temperature must be finite.");

report = struct;
report.sourceFile = filePath;
report.selectedSourceRows = selectedSourceRows;
report.excludedSourceRows = setdiff((1:sampleCount)', selectedSourceRows, 'stable');
report.removedDuplicateSourceRows = removedDuplicateSourceRows;
report.keptSourceRows = keptSourceRows;
report.currentNormalisation = "I = -meas.Current (applied once)";
report.interpolationApplied = false;
end

function isDuplicate = rowsAreExactDuplicates(meas, requiredFields, rowA, rowB)
isDuplicate = true;
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields(fieldIndex);
    if fieldName == "TimeStamp"
        valuesMatch = string(meas.(fieldName)(rowA)) == ...
            string(meas.(fieldName)(rowB));
    else
        valuesMatch = isequaln(meas.(fieldName)(rowA), meas.(fieldName)(rowB));
    end
    isDuplicate = isDuplicate && valuesMatch;
end
end
