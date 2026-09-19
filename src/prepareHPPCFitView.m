function fitView = prepareHPPCFitView(meas)
%PREPAREHPPCFITVIEW Build the strict-time HPPC fit view.
%   FITVIEW = PREPAREHPPCFITVIEW(MEAS) takes the raw HPPC meas struct and
%   returns a derived view that still tracks source-row provenance. Drops
%   the second row of exact zero-dt duplicates, excludes both rows of a
%   nonidentical zero-dt pair (can't tell which one is right), and
%   otherwise keeps every raw row and its actual numeric interval - no
%   interpolation, no fudging the timestamps.
%
%   Current gets normalised once here:
%     FITVIEW.I = -MEAS.Current
%   so positive means discharge. This function just prepares the view,
%   no state updates, no parameter fitting, no picking "representative"
%   rows.

arguments
    meas (1, 1) struct
end

requiredFields = [ ...
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

assert(~isempty(fieldnames(meas)), "prepareHPPCFitView:InvalidMeas", ...
    "The raw HPPC meas struct must not be empty.");
missingFields = setdiff(requiredFields, string(fieldnames(meas)));
assert(isempty(missingFields), "prepareHPPCFitView:MissingFields", ...
    "The raw HPPC meas struct is missing: %s.", join(missingFields, ", "));

rawSampleCount = numel(meas.Time);
assert(rawSampleCount > 0, "prepareHPPCFitView:EmptyInput", ...
    "The raw HPPC meas struct must contain at least one sample.");
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields(fieldIndex);
    assert(isvector(meas.(fieldName)) && ...
        numel(meas.(fieldName)) == rawSampleCount, ...
        "prepareHPPCFitView:InconsistentFieldLength", ...
        "Field %s must be a vector with %d samples.", fieldName, rawSampleCount);
end

numericFields = requiredFields(requiredFields ~= "TimeStamp");
for fieldIndex = 1:numel(numericFields)
    fieldName = numericFields(fieldIndex);
    assert(isnumeric(meas.(fieldName)) && all(isfinite(meas.(fieldName))), ...
        "prepareHPPCFitView:NonFiniteSignal", ...
        "Numeric field %s must contain only finite values.", fieldName);
end

rawTime_s = meas.Time(:);
dt_s = diff(rawTime_s);
assert(~any(dt_s < 0), "prepareHPPCFitView:DecreasingTime", ...
    "Raw HPPC numeric Time must not decrease.");

zeroDtSecondRows = find(dt_s == 0) + 1;
exactDuplicateRemovedSourceRows = zeros(0, 1);
nonidenticalExcludedSourceRows = zeros(0, 1);
for secondRow = zeroDtSecondRows.'
    firstRow = secondRow - 1;
    if rowsAreExactDuplicates(meas, requiredFields, firstRow, secondRow)
        exactDuplicateRemovedSourceRows(end + 1, 1) = secondRow; %#ok<SAGROW>
    else
        nonidenticalExcludedSourceRows(end + 1, 1) = firstRow; %#ok<SAGROW>
        nonidenticalExcludedSourceRows(end + 1, 1) = secondRow; %#ok<SAGROW>
    end
end

assert(numel(unique(nonidenticalExcludedSourceRows)) == ...
    numel(nonidenticalExcludedSourceRows), ...
    "prepareHPPCFitView:AmbiguousZeroDtRun", ...
    "A raw row belongs to more than one nonidentical zero-dt pair.");
assert(isempty(intersect(exactDuplicateRemovedSourceRows, ...
    nonidenticalExcludedSourceRows)), ...
    "prepareHPPCFitView:AmbiguousZeroDtRun", ...
    "An exact-duplicate removal overlaps a nonidentical zero-dt pair.");

removeMask = false(rawSampleCount, 1);
removeMask(exactDuplicateRemovedSourceRows) = true;
removeMask(nonidenticalExcludedSourceRows) = true;
retainedSourceRows = find(~removeMask);
derivedTime_s = rawTime_s(retainedSourceRows);
assert(all(diff(derivedTime_s) > 0), ...
    "prepareHPPCFitView:InvalidDerivedTime", ...
    "The derived HPPC fit view must have strictly increasing numeric time.");

fitView = struct;
fitView.t = derivedTime_s;
fitView.I = -meas.Current(retainedSourceRows);
fitView.V = meas.Voltage(retainedSourceRows);
fitView.rawCurrent = meas.Current(retainedSourceRows);
fitView.Ah = meas.Ah(retainedSourceRows);
fitView.Wh = meas.Wh(retainedSourceRows);
fitView.Power = meas.Power(retainedSourceRows);
fitView.T = meas.Battery_Temp_degC(retainedSourceRows);
fitView.chamberTemperature_degC = meas.Chamber_Temp_degC(retainedSourceRows);
fitView.sourceRow = retainedSourceRows;

fitView.report = struct;
fitView.report.rawSampleCount = rawSampleCount;
fitView.report.retainedSourceRows = retainedSourceRows;
fitView.report.exactDuplicateRemovedSourceRows = ...
    exactDuplicateRemovedSourceRows;
fitView.report.nonidenticalExcludedSourceRows = ...
    nonidenticalExcludedSourceRows;
fitView.report.currentNormalisation = "I = -meas.Current (applied once)";
fitView.report.interpolationApplied = false;
fitView.report.timeRepairApplied = false;
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
